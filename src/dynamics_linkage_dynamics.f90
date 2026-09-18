module dynamics_linkage_dynamics
	!! Provides a maximal-coordinate dynamics adapter for serial and parallel
	!! linkages. Link mass properties and joint attachment frames are converted
	!! into the rigid bodies and holonomic constraints required by the
	!! variational integrator.
	use iso_fortran_env, only : int32, real64
	use dynamics_error_handling, only : DYN_ARRAY_SIZE_ERROR, &
		DYN_INVALID_INPUT_ERROR
	use dynamics_kinematics, only : transform_inverse
	use dynamics_linkage, only : link, binary_link, serial_linkage
	use dynamics_parallel_linkage, only : kinematic_mechanism, planar_linkage
	use dynamics_joints, only : joint, FIXED_JOINT, REVOLUTE_JOINT, &
		PRISMATIC_JOINT, CYLINDRICAL_JOINT, UNIVERSAL_JOINT, SPHERICAL_JOINT
	use dynamics_quaternions, only : quaternion
	use dynamics_rigid_bodies, only : rigid_body
	use dynamics_variational_integrators, only : variational_integrator, &
		variational_state
	use linalg, only : identity
	implicit none
	private

	public :: linkage_dynamic_model
	public :: linkage_prescribed_motion
	public :: joint_reaction

	abstract interface
		function linkage_prescribed_motion(t) result(rst)
			!! Computes a prescribed angular displacement as a function of time.
			import :: real64
			real(real64), intent(in) :: t
				!! The simulation time.
			real(real64) :: rst
				!! The prescribed absolute angle, in radians.
		end function
	end interface

	type dynamic_joint
		!! Stores one joint in maximal-coordinate body indexing.
		integer(int32) :: joint_type = REVOLUTE_JOINT
		integer(int32) :: parent_body = 0
		integer(int32) :: child_body = 0
		real(real64), dimension(4,4) :: parent_frame
		real(real64), dimension(4,4) :: child_frame
	end type

	type joint_reaction
		!! Defines the constraint reaction exerted by a joint on its child link.
		!! Both vectors are expressed in the world coordinate frame. The reaction
		!! exerted on the parent link is equal and opposite.
		real(real64), dimension(3) :: force = 0.0d0
			!! The joint reaction force.
		real(real64), dimension(3) :: moment = 0.0d0
			!! The joint reaction moment about the joint center.
	end type

	type linkage_dynamic_model
		!! Defines the variational-integrator representation of a linkage.
		!! Body index zero in a joint descriptor denotes the fixed ground link.
		type(rigid_body), allocatable, private, dimension(:) :: m_bodies
		type(dynamic_joint), allocatable, private, dimension(:) :: m_joints
		type(variational_state), private :: m_initial_state
		logical, private :: m_planar = .false.
		integer(int32), private :: m_constraint_count = 0
	contains
		procedure, public :: get_body_count => ldm_get_body_count
		procedure, public :: get_joint_count => ldm_get_joint_count
		procedure, public :: get_constraint_count => ldm_get_constraint_count
		procedure, public :: get_initial_state => ldm_get_initial_state
		procedure, public :: constraint_residual => ldm_constraint_residual
		procedure, public :: get_joint_reactions => ldm_get_joint_reactions
		procedure, public :: solve => ldm_solve
	end type

	interface linkage_dynamic_model
		module procedure :: ldm_from_serial
		module procedure :: ldm_from_parallel
	end interface

	type linkage_solve_context
		class(linkage_dynamic_model), pointer :: model => null()
		real(real64), dimension(3) :: gravity = 0.0d0
		real(real64), allocatable, dimension(:,:) :: body_force
		real(real64), allocatable, dimension(:,:) :: body_torque
		integer(int32) :: prescribed_body = 0
		procedure(linkage_prescribed_motion), pointer, nopass :: &
			prescribed_motion => null()
	end type

contains
! ------------------------------------------------------------------------------
function ldm_from_serial(mechanism, q) result(rst)
	!! Constructs a dynamic model from a serial linkage and a compatible set of
	!! joint variables. Every serial link is treated as a moving rigid body;
	!! the proximal joint of the first link is attached to ground.
	type(serial_linkage), intent(in) :: mechanism
		!! The serial linkage to convert.
	real(real64), intent(in), dimension(:) :: q
		!! One joint variable for each link.
	type(linkage_dynamic_model) :: rst
		!! The resulting dynamic model.

	integer(int32) :: i, n
	real(real64), dimension(4,4) :: body_transform, motion
	class(binary_link), pointer :: current_link

	n = mechanism%get_link_count()
	if (size(q) /= n) error stop DYN_ARRAY_SIZE_ERROR
	allocate(rst%m_bodies(n), rst%m_joints(n))
	call allocate_state(rst%m_initial_state, n)
	body_transform = identity(4)
	do i = 1, n
		current_link => mechanism%get_link(i)
		call copy_mass_properties(rst%m_bodies(i), current_link)
		rst%m_joints(i)%joint_type = current_link%joint_type
		rst%m_joints(i)%parent_body = i - 1
		rst%m_joints(i)%child_body = i
		rst%m_joints(i)%parent_frame = identity(4)
		rst%m_joints(i)%child_frame = current_link%get_joint_frame(1)
		motion = serial_joint_motion(current_link%joint_type, q(i))
		body_transform = matmul(body_transform, matmul(motion, &
			transform_inverse(rst%m_joints(i)%child_frame)))
		call set_body_state(rst%m_initial_state, i, current_link, body_transform)
	end do
	rst%m_constraint_count = sum_joint_constraints(rst%m_joints, .false.)
end function

! ------------------------------------------------------------------------------
function ldm_from_parallel(mechanism, q) result(rst)
	!! Constructs a dynamic model from a parallel or planar linkage. The
	!! mechanism's base link remains fixed and all other links become moving
	!! maximal-coordinate rigid bodies.
	class(kinematic_mechanism), intent(in) :: mechanism
		!! The parallel linkage to convert.
	real(real64), intent(in), dimension(:) :: q
		!! The complete, constraint-compatible joint-variable array.
	type(linkage_dynamic_model) :: rst
		!! The resulting dynamic model.

	integer(int32) :: base, body_index, i, nbody, nlink, njoint
	integer(int32), allocatable, dimension(:) :: body_map
	real(real64), dimension(4,4) :: body_transform
	type(joint) :: current_joint
	class(link), pointer :: current_link

	if (size(q) /= mechanism%get_variable_count()) &
		error stop DYN_ARRAY_SIZE_ERROR
	nlink = mechanism%get_link_count()
	njoint = mechanism%get_joint_count()
	base = mechanism%get_base_link()
	nbody = nlink - 1
	if (nbody < 1) error stop DYN_INVALID_INPUT_ERROR
	allocate(rst%m_bodies(nbody), rst%m_joints(njoint), body_map(nlink))
	call allocate_state(rst%m_initial_state, nbody)
	body_map = 0
	body_index = 0
	do i = 1, nlink
		if (i == base) cycle
		body_index = body_index + 1
		body_map(i) = body_index
		current_link => mechanism%get_link(i)
		call copy_mass_properties(rst%m_bodies(body_index), current_link)
		body_transform = mechanism%body_transform(i, q)
		call set_body_state(rst%m_initial_state, body_index, current_link, &
			body_transform)
	end do

	do i = 1, njoint
		current_joint = mechanism%get_joint(i)
		rst%m_joints(i)%joint_type = current_joint%joint_type
		rst%m_joints(i)%parent_body = body_map(current_joint%parent_link)
		rst%m_joints(i)%child_body = body_map(current_joint%child_link)
		current_link => mechanism%get_link(current_joint%parent_link)
		rst%m_joints(i)%parent_frame = &
			current_link%get_joint_frame(current_joint%parent_frame)
		current_link => mechanism%get_link(current_joint%child_link)
		rst%m_joints(i)%child_frame = &
			current_link%get_joint_frame(current_joint%child_frame)
	end do
	select type (mechanism)
	type is (planar_linkage)
		rst%m_planar = .true.
	class default
		rst%m_planar = .false.
	end select
	rst%m_constraint_count = sum_joint_constraints(rst%m_joints, &
		rst%m_planar)
	if (rst%m_planar) rst%m_constraint_count = &
		rst%m_constraint_count + 3 * nbody
end function

! ------------------------------------------------------------------------------
pure function ldm_get_body_count(this) result(rst)
	class(linkage_dynamic_model), intent(in) :: this
	integer(int32) :: rst
	rst = size(this%m_bodies)
end function

! ------------------------------------------------------------------------------
pure function ldm_get_joint_count(this) result(rst)
	!! Gets the number of joints represented by the dynamic model.
	class(linkage_dynamic_model), intent(in) :: this
		!! The dynamic linkage model.
	integer(int32) :: rst
		!! The joint count.

	rst = size(this%m_joints)
end function

! ------------------------------------------------------------------------------
pure function ldm_get_constraint_count(this) result(rst)
	class(linkage_dynamic_model), intent(in) :: this
	integer(int32) :: rst
	rst = this%m_constraint_count
end function

! ------------------------------------------------------------------------------
function ldm_get_initial_state(this) result(rst)
	!! Gets a copy of the constraint-compatible state used to construct the
	!! dynamic model.
	class(linkage_dynamic_model), intent(in) :: this
	type(variational_state) :: rst
	rst = this%m_initial_state
end function

! ------------------------------------------------------------------------------
function ldm_constraint_residual(this, state) result(rst)
	!! Evaluates every joint and planar constraint for a supplied state.
	class(linkage_dynamic_model), intent(in) :: this
	type(variational_state), intent(in) :: state
	real(real64), allocatable, dimension(:) :: rst
	integer(int32) :: i, index

	allocate(rst(this%m_constraint_count))
	index = 1
	if (this%m_planar) then
		do i = 1, this%get_body_count()
			call append_planar_body_constraints(this, state, i, rst, index)
		end do
	end if
	do i = 1, size(this%m_joints)
		call append_joint_constraints(this, state, this%m_joints(i), rst, index)
	end do
end function

! ------------------------------------------------------------------------------
function ldm_get_joint_reactions(this, state, multipliers) result(rst)
	!! Converts one time step's constraint multipliers into the force and moment
	!! exerted by every joint on its child link. Multipliers associated with the
	!! planar-body constraints or a prescribed-motion constraint are excluded.
	class(linkage_dynamic_model), intent(in) :: this
		!! The dynamic linkage model.
	type(variational_state), intent(in) :: state
		!! The state corresponding to the supplied multipliers.
	real(real64), intent(in), dimension(:) :: multipliers
		!! The constraint multiplier vector returned by the integrator.
	type(joint_reaction), allocatable, dimension(:) :: rst
		!! One world-frame reaction wrench for each joint, in mechanism order.

	integer(int32) :: i, index, nconstraint
	real(real64), dimension(3) :: local_force, local_moment
	real(real64), dimension(4,4) :: parent

	if (size(multipliers) < this%m_constraint_count) &
		error stop DYN_ARRAY_SIZE_ERROR
	allocate(rst(size(this%m_joints)))
	index = 1
	if (this%m_planar) index = 3 * this%get_body_count() + 1
	do i = 1, size(this%m_joints)
		nconstraint = joint_constraint_count(this%m_joints(i)%joint_type, &
			this%m_planar)
		if (this%m_planar) then
			call planar_joint_reaction(this, state, this%m_joints(i), &
				multipliers(index:index+nconstraint-1), rst(i))
		else
			call spatial_joint_reaction(this%m_joints(i)%joint_type, &
				multipliers(index:index+nconstraint-1), local_force, local_moment)
			parent = joint_world_transform(this, state, &
				this%m_joints(i)%parent_body, this%m_joints(i)%parent_frame)
			rst(i)%force = matmul(parent(1:3,1:3), local_force)
			rst(i)%moment = matmul(parent(1:3,1:3), local_moment)
		end if
		index = index + nconstraint
	end do
end function

! ------------------------------------------------------------------------------
subroutine planar_joint_reaction(model, state, descriptor, multipliers, reaction)
	!! Maps planar joint multipliers to a world-frame reaction wrench.
	class(linkage_dynamic_model), intent(in) :: model
	type(variational_state), intent(in) :: state
	type(dynamic_joint), intent(in) :: descriptor
	real(real64), intent(in), dimension(:) :: multipliers
	type(joint_reaction), intent(out) :: reaction

	real(real64), dimension(4,4) :: parent
	real(real64), dimension(3) :: axis, normal

	reaction%force = 0.0d0
	reaction%moment = 0.0d0
	select case (descriptor%joint_type)
	case (FIXED_JOINT)
		reaction%force(1:2) = multipliers(1:2)
		reaction%moment(3) = multipliers(3)
	case (REVOLUTE_JOINT)
		reaction%force(1:2) = multipliers(1:2)
	case (PRISMATIC_JOINT)
		parent = joint_world_transform(model, state, descriptor%parent_body, &
			descriptor%parent_frame)
		axis = parent(1:3,3)
		normal = [-axis(2), axis(1), 0.0d0]
		reaction%force = multipliers(1) * normal
		reaction%moment(3) = multipliers(2)
	case default
		error stop DYN_INVALID_INPUT_ERROR
	end select
end subroutine

! ------------------------------------------------------------------------------
pure subroutine spatial_joint_reaction(joint_type, multipliers, force, moment)
	!! Maps spatial joint multipliers into a wrench expressed in the parent
	!! joint frame. The caller rotates the wrench into world coordinates.
	integer(int32), intent(in) :: joint_type
	real(real64), intent(in), dimension(:) :: multipliers
	real(real64), intent(out), dimension(3) :: force, moment

	force = 0.0d0
	moment = 0.0d0
	select case (joint_type)
	case (FIXED_JOINT)
		force = multipliers(1:3)
		moment = multipliers(4:6)
	case (REVOLUTE_JOINT)
		force = multipliers(1:3)
		moment(1) = multipliers(5)
		moment(2) = -multipliers(4)
	case (PRISMATIC_JOINT)
		force(1:2) = multipliers(1:2)
		moment = multipliers(3:5)
	case (CYLINDRICAL_JOINT)
		force(1:2) = multipliers(1:2)
		moment(1) = multipliers(4)
		moment(2) = -multipliers(3)
	case (UNIVERSAL_JOINT)
		force = multipliers(1:3)
		moment(3) = -multipliers(4)
	case (SPHERICAL_JOINT)
		force = multipliers(1:3)
	case default
		error stop DYN_INVALID_INPUT_ERROR
	end select
end subroutine

! ------------------------------------------------------------------------------
function ldm_solve(this, integrator, dt, ntime, initial_state, gravity, &
	body_force, body_torque, prescribed_body, prescribed_motion, multipliers) &
	result(rst)
	!! Integrates the linkage dynamics under an optional uniform world-frame
	!! gravitational acceleration and optional constant body loads.
	class(linkage_dynamic_model), intent(in), target :: this
	type(variational_integrator), intent(in) :: integrator
	real(real64), intent(in) :: dt
	integer(int32), intent(in) :: ntime
	type(variational_state), intent(in), optional :: initial_state
	real(real64), intent(in), optional, dimension(3) :: gravity
	real(real64), intent(in), optional, dimension(:,:) :: body_force
		!! Constant 3-by-nbody world-frame force array.
	real(real64), intent(in), optional, dimension(:,:) :: body_torque
		!! Constant 3-by-nbody body-frame torque array.
	integer(int32), intent(in), optional :: prescribed_body
		!! The moving body whose absolute planar angle is prescribed.
	procedure(linkage_prescribed_motion), optional :: prescribed_motion
		!! The prescribed absolute planar angle as a function of time.
	real(real64), allocatable, intent(out), optional, dimension(:,:) :: multipliers
		!! Constraint multipliers for each completed time step. When a motion is
		!! prescribed, the last row is the required motor torque.
	type(variational_state), allocatable, dimension(:) :: rst

	type(linkage_solve_context) :: context
	type(variational_state) :: state
	integer(int32) :: constraint_count

	state = this%m_initial_state
	if (present(initial_state)) state = initial_state
	context%model => this
	if (present(gravity)) context%gravity = gravity
	allocate(context%body_force(3,this%get_body_count()), source = 0.0d0)
	allocate(context%body_torque(3,this%get_body_count()), source = 0.0d0)
	if (present(body_force)) then
		if (any(shape(body_force) /= [3, this%get_body_count()])) &
			error stop DYN_ARRAY_SIZE_ERROR
		context%body_force = body_force
	end if
	if (present(body_torque)) then
		if (any(shape(body_torque) /= [3, this%get_body_count()])) &
			error stop DYN_ARRAY_SIZE_ERROR
		context%body_torque = body_torque
	end if
	if (present(prescribed_body) .neqv. present(prescribed_motion)) &
		error stop DYN_INVALID_INPUT_ERROR
	constraint_count = this%m_constraint_count
	if (present(prescribed_body)) then
		if (.not.this%m_planar .or. prescribed_body < 1 .or. &
			prescribed_body > this%get_body_count()) &
			error stop DYN_INVALID_INPUT_ERROR
		context%prescribed_body = prescribed_body
		context%prescribed_motion => prescribed_motion
		constraint_count = constraint_count + 1
	end if
	rst = integrator%solve(this%m_bodies, state, dt, ntime, &
		constraint_count = constraint_count, &
		constraint = linkage_constraints, &
		force_function = linkage_gravity, multipliers = multipliers, &
		args = context)
end function

! ------------------------------------------------------------------------------
subroutine linkage_constraints(state, value, args)
	type(variational_state), intent(in) :: state
	real(real64), intent(out), dimension(:) :: value
	class(*), intent(inout), optional :: args

	select type (context => args)
	type is (linkage_solve_context)
		value(1:context%model%m_constraint_count) = &
			context%model%constraint_residual(state)
		if (associated(context%prescribed_motion)) then
			value(size(value)) = planar_body_angle(state, &
				context%prescribed_body) - &
				context%prescribed_motion(state%time)
		end if
	class default
		error stop DYN_INVALID_INPUT_ERROR
	end select
end subroutine

! ------------------------------------------------------------------------------
function planar_body_angle(state, body) result(rst)
	!! Gets the absolute angle of a planar body about the world z-axis.
	type(variational_state), intent(in) :: state
	integer(int32), intent(in) :: body
	real(real64) :: rst
	real(real64), dimension(3,3) :: rotation

	rotation = state%orientation(body)%to_matrix()
	rst = atan2(rotation(2,1), rotation(1,1))
end function

! ------------------------------------------------------------------------------
subroutine linkage_gravity(t, state, force, torque, args)
	real(real64), intent(in) :: t
	type(variational_state), intent(in) :: state
	real(real64), intent(out), dimension(:,:) :: force, torque
	class(*), intent(inout), optional :: args
	integer(int32) :: i

	force = 0.0d0
	torque = 0.0d0
	select type (context => args)
	type is (linkage_solve_context)
		do i = 1, context%model%get_body_count()
			force(:,i) = context%model%m_bodies(i)%mass * context%gravity
		end do
		force = force + context%body_force
		torque = torque + context%body_torque
	class default
		error stop DYN_INVALID_INPUT_ERROR
	end select
end subroutine

! ------------------------------------------------------------------------------
subroutine append_planar_body_constraints(model, state, body, value, index)
	class(linkage_dynamic_model), intent(in) :: model
	type(variational_state), intent(in) :: state
	integer(int32), intent(in) :: body
	real(real64), intent(inout), dimension(:) :: value
	integer(int32), intent(inout) :: index
	real(real64), dimension(3,3) :: rotation
	real(real64), dimension(3) :: origin

	call body_pose(model, state, body, rotation, origin)
	value(index) = origin(3)
	value(index+1) = rotation(3,1)
	value(index+2) = rotation(3,2)
	index = index + 3
end subroutine

! ------------------------------------------------------------------------------
subroutine append_joint_constraints(model, state, descriptor, value, index)
	class(linkage_dynamic_model), intent(in) :: model
	type(variational_state), intent(in) :: state
	type(dynamic_joint), intent(in) :: descriptor
	real(real64), intent(inout), dimension(:) :: value
	integer(int32), intent(inout) :: index
	real(real64), dimension(4,4) :: parent, child, relative
	real(real64), dimension(3) :: displacement, axis, normal
	real(real64) :: angle

	parent = joint_world_transform(model, state, descriptor%parent_body, &
		descriptor%parent_frame)
	child = joint_world_transform(model, state, descriptor%child_body, &
		descriptor%child_frame)
	relative = matmul(transform_inverse(parent), child)
	if (model%m_planar) then
		displacement = child(1:3,4) - parent(1:3,4)
		angle = atan2(relative(2,1), relative(1,1))
		select case (descriptor%joint_type)
		case (FIXED_JOINT)
			value(index:index+1) = displacement(1:2)
			value(index+2) = angle
			index = index + 3
		case (REVOLUTE_JOINT)
			value(index:index+1) = displacement(1:2)
			index = index + 2
		case (PRISMATIC_JOINT)
			axis = parent(1:3,3)
			normal = [-axis(2), axis(1), 0.0d0]
			value(index) = dot_product(displacement, normal)
			value(index+1) = angle
			index = index + 2
		case default
			error stop DYN_INVALID_INPUT_ERROR
		end select
	else
		call append_spatial_joint_constraints(descriptor%joint_type, relative, &
			value, index)
	end if
end subroutine

! ------------------------------------------------------------------------------
subroutine append_spatial_joint_constraints(joint_type, relative, value, index)
	integer(int32), intent(in) :: joint_type
	real(real64), intent(in), dimension(4,4) :: relative
	real(real64), intent(inout), dimension(:) :: value
	integer(int32), intent(inout) :: index
	real(real64), dimension(3) :: rotation_error

	rotation_error = 0.5d0 * [relative(3,2) - relative(2,3), &
		relative(1,3) - relative(3,1), relative(2,1) - relative(1,2)]
	select case (joint_type)
	case (FIXED_JOINT)
		value(index:index+2) = relative(1:3,4)
		value(index+3:index+5) = rotation_error
		index = index + 6
	case (REVOLUTE_JOINT)
		value(index:index+2) = relative(1:3,4)
		value(index+3:index+4) = relative(3,1:2)
		index = index + 5
	case (PRISMATIC_JOINT)
		value(index:index+1) = relative(1:2,4)
		value(index+2:index+4) = rotation_error
		index = index + 5
	case (CYLINDRICAL_JOINT)
		value(index:index+1) = relative(1:2,4)
		value(index+2:index+3) = relative(3,1:2)
		index = index + 4
	case (UNIVERSAL_JOINT)
		value(index:index+2) = relative(1:3,4)
		value(index+3) = relative(1,2)
		index = index + 4
	case (SPHERICAL_JOINT)
		value(index:index+2) = relative(1:3,4)
		index = index + 3
	case default
		error stop DYN_INVALID_INPUT_ERROR
	end select
end subroutine

! ------------------------------------------------------------------------------
function joint_world_transform(model, state, body, frame) result(rst)
	class(linkage_dynamic_model), intent(in) :: model
	type(variational_state), intent(in) :: state
	integer(int32), intent(in) :: body
	real(real64), intent(in), dimension(4,4) :: frame
	real(real64), dimension(4,4) :: rst, transform
	real(real64), dimension(3,3) :: rotation
	real(real64), dimension(3) :: origin

	if (body == 0) then
		rst = frame
	else
		call body_pose(model, state, body, rotation, origin)
		transform = identity(4)
		transform(1:3,1:3) = rotation
		transform(1:3,4) = origin
		rst = matmul(transform, frame)
	end if
end function

! ------------------------------------------------------------------------------
subroutine body_pose(model, state, body, rotation, origin)
	class(linkage_dynamic_model), intent(in) :: model
	type(variational_state), intent(in) :: state
	integer(int32), intent(in) :: body
	real(real64), intent(out), dimension(3,3) :: rotation
	real(real64), intent(out), dimension(3) :: origin

	rotation = state%orientation(body)%to_matrix()
	origin = state%position(:,body) - &
		matmul(rotation, model%m_bodies(body)%cg)
end subroutine

! ------------------------------------------------------------------------------
subroutine allocate_state(state, nbody)
	type(variational_state), intent(out) :: state
	integer(int32), intent(in) :: nbody
	integer(int32) :: i

	allocate(state%position(3,nbody), state%orientation(nbody), &
		state%velocity(3,nbody), state%angular_velocity(3,nbody))
	state%position = 0.0d0
	state%velocity = 0.0d0
	state%angular_velocity = 0.0d0
	state%time = 0.0d0
	do i = 1, nbody
		state%orientation(i) = quaternion([1.0d0, 0.0d0, 0.0d0, 0.0d0])
	end do
end subroutine

! ------------------------------------------------------------------------------
subroutine set_body_state(state, body, properties, transform)
	type(variational_state), intent(inout) :: state
	integer(int32), intent(in) :: body
	class(link), intent(in) :: properties
	real(real64), intent(in), dimension(4,4) :: transform

	state%orientation(body) = quaternion(transform(1:3,1:3))
	state%position(:,body) = transform(1:3,4) + &
		matmul(transform(1:3,1:3), properties%cg)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine copy_mass_properties(body, source)
	!! Copies the inherited rigid-body properties from a linkage link.
	type(rigid_body), intent(out) :: body
	class(link), intent(in) :: source

	body%mass = source%mass
	body%cg = source%cg
	body%inertia = source%inertia
end subroutine

! ------------------------------------------------------------------------------
pure function serial_joint_motion(joint_type, q) result(rst)
	integer(int32), intent(in) :: joint_type
	real(real64), intent(in) :: q
	real(real64), dimension(4,4) :: rst
	real(real64) :: c, s

	rst = identity(4)
	select case (joint_type)
	case (REVOLUTE_JOINT)
		c = cos(q)
		s = sin(q)
		rst(1,1) = c
		rst(2,1) = s
		rst(1,2) = -s
		rst(2,2) = c
	case (PRISMATIC_JOINT)
		rst(3,4) = q
	case default
		error stop DYN_INVALID_INPUT_ERROR
	end select
end function

! ------------------------------------------------------------------------------
pure function sum_joint_constraints(joints, planar) result(rst)
	type(dynamic_joint), intent(in), dimension(:) :: joints
	logical, intent(in) :: planar
	integer(int32) :: rst
	integer(int32) :: i

	rst = 0
	do i = 1, size(joints)
		rst = rst + joint_constraint_count(joints(i)%joint_type, planar)
	end do
end function

! ------------------------------------------------------------------------------
pure function joint_constraint_count(joint_type, planar) result(rst)
	!! Gets the number of scalar multipliers associated with one joint.
	integer(int32), intent(in) :: joint_type
	logical, intent(in) :: planar
	integer(int32) :: rst

	if (planar) then
		select case (joint_type)
		case (FIXED_JOINT)
			rst = 3
		case (REVOLUTE_JOINT, PRISMATIC_JOINT)
			rst = 2
		case default
			error stop DYN_INVALID_INPUT_ERROR
		end select
	else
		select case (joint_type)
		case (FIXED_JOINT)
			rst = 6
		case (REVOLUTE_JOINT, PRISMATIC_JOINT)
			rst = 5
		case (CYLINDRICAL_JOINT, UNIVERSAL_JOINT)
			rst = 4
		case (SPHERICAL_JOINT)
			rst = 3
		case default
			error stop DYN_INVALID_INPUT_ERROR
		end select
	end if
end function

end module

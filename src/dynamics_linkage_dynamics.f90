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
    use dynamics_helper, only : cross_product
    use dynamics_rigid_bodies, only : rigid_body
    use dynamics_variational_integrators
    use linalg, only : identity
    implicit none
    private
    public :: linkage_dynamic_model
    public :: linkage_prescribed_motion
    public :: joint_reaction
    public :: dynamic_joint
    public :: axial_force_element
    public :: linear_spring
    public :: linear_damper
    public :: torsional_force_element
    public :: torsional_spring
    public :: torsional_damper
    public :: axial_element_result
    public :: torsional_element_result

    abstract interface
        function linkage_prescribed_motion(t, args) result(rst)
            !! Computes a prescribed angular displacement as a function of time.
            import :: real64
            real(real64), intent(in) :: t
                !! The simulation time.
			class(*), intent(inout), optional :: args
				!! A mechanism for passing information in/out of this routine.
            real(real64) :: rst
                !! The prescribed absolute angle, in radians.
        end function
    end interface

    type dynamic_joint
        !! Describes one joint constraint in maximal-coordinate body indexing.
        !! Body index zero denotes the fixed ground link. Each frame is a
        !! homogeneous transform from its associated body frame to the joint
        !! frame; for a ground body, the frame is expressed directly in world
        !! coordinates. The parent and child joint frames coincide when the
        !! joint constraint is satisfied.
        integer(int32) :: joint_type = REVOLUTE_JOINT
            !! The joint kind, using one of the joint-type constants from
            !! dynamics_joints.
        integer(int32) :: parent_body = 0
            !! The one-based maximal-coordinate body index on the parent side,
            !! or zero for the fixed ground link.
        integer(int32) :: child_body = 0
            !! The one-based maximal-coordinate body index on the child side,
            !! or zero for the fixed ground link.
        real(real64), dimension(4,4) :: parent_frame
            !! Homogeneous transform locating and orienting the joint frame in
            !! the parent body frame, or in world coordinates when parent_body is
            !! zero.
        real(real64), dimension(4,4) :: child_frame
            !! Homogeneous transform locating and orienting the joint frame in
            !! the child body frame, or in world coordinates when child_body is
            !! zero.
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

    type axial_element_result
        !! Defines the instantaneous state and scalar force of an axial element.
        real(real64) :: length = 0.0d0
            !! The current distance between attachment points.
        real(real64) :: length_rate = 0.0d0
            !! The relative attachment velocity along the element axis.
        real(real64) :: force = 0.0d0
            !! The signed force acting on body 1 toward body 2.
    end type

    type torsional_element_result
        !! Defines the instantaneous twist state and torque of a torsional element.
        real(real64) :: angle = 0.0d0
            !! The signed joint angle, in radians.
        real(real64) :: angle_rate = 0.0d0
            !! The relative twist rate about the revolute-joint axis.
        real(real64) :: torque = 0.0d0
            !! The signed torque applied to the child body about the joint axis.
    end type

    type axial_force_element
        !! Defines an extensible force element between two attachment points.
        !! The element force acts along the line from point 1 to point 2.
        !! Body index zero denotes a fixed ground point whose coordinates are
        !! expressed in world coordinates. Subtypes override evaluate_force to
        !! compute a force from the current length and length rate.
        integer(int32) :: body_1 = 0
            !! The one-based body index for the first attachment point, or zero
            !! for a ground point.
        integer(int32) :: body_2 = 0
            !! The one-based body index for the second attachment point, or zero
            !! for a ground point.
        real(real64), dimension(3) :: point_1 = 0.0d0
            !! The first attachment point in body_1 coordinates, or in world
            !! coordinates when body_1 is zero.
        real(real64), dimension(3) :: point_2 = 0.0d0
            !! The second attachment point in body_2 coordinates, or in world
            !! coordinates when body_2 is zero.
    contains
        procedure, public :: evaluate_force => axial_zero_force
            !! Returns the signed axial force for a given element length and
            !! length rate. The base implementation returns zero.
    end type

    type, extends(axial_force_element) :: linear_spring
        !! Defines a linear axial spring supporting tension and compression.
        real(real64) :: stiffness = 0.0d0
            !! The force per unit extension.
        real(real64) :: free_length = 0.0d0
            !! The zero-force element length, including preload definition.
    contains
        procedure, public :: evaluate_force => linear_spring_force
    end type

    type, extends(axial_force_element) :: linear_damper
        !! Defines a linear viscous damper acting only along the element axis.
        real(real64) :: damping = 0.0d0
            !! The force per unit axial relative velocity.
    contains
        procedure, public :: evaluate_force => linear_damper_force
    end type

    type torsional_force_element
        !! Defines an extensible torsional element bound to a revolute joint.
        integer(int32) :: joint_index = 0
            !! The one-based revolute-joint index whose axis the element uses.
    contains
        procedure, public :: evaluate_torque => torsional_zero_torque
    end type

    type, extends(torsional_force_element) :: torsional_spring
        !! Defines a linear torsional spring about a revolute-joint axis.
        real(real64) :: stiffness = 0.0d0
            !! The torque per unit angular displacement.
        real(real64) :: free_angle = 0.0d0
            !! The zero-torque relative joint angle, in radians.
    contains
        procedure, public :: evaluate_torque => torsional_spring_torque
    end type

    type, extends(torsional_force_element) :: torsional_damper
        !! Defines a linear damper opposing only joint-axis twist rate.
        real(real64) :: damping = 0.0d0
            !! The torque per unit relative angular velocity.
    contains
        procedure, public :: evaluate_torque => torsional_damper_torque
    end type

    type axial_element_container
        !! Holds one polymorphic axial force element in the dynamic model.
        class(axial_force_element), allocatable :: item
            !! The allocated axial element stored by the container.
    end type

    type torsional_element_container
        !! Holds one polymorphic torsional force element in the dynamic model.
        class(torsional_force_element), allocatable :: item
            !! The allocated torsional element stored by the container.
    end type

    type linkage_dynamic_model
        !! Defines the variational-integrator representation of a linkage.
        !! Body index zero in a joint descriptor denotes the fixed ground link.
        type(rigid_body), allocatable, private, dimension(:) :: m_bodies
        type(dynamic_joint), allocatable, private, dimension(:) :: m_joints
        type(axial_element_container), allocatable, private, dimension(:) :: &
            m_axial_elements
        type(torsional_element_container), allocatable, private, dimension(:) :: &
            m_torsional_elements
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
        procedure, public :: add_axial_element => ldm_add_axial_element
        procedure, public :: add_linear_spring => ldm_add_linear_spring
        procedure, public :: add_linear_damper => ldm_add_linear_damper
        procedure, public :: add_torsional_element => ldm_add_torsional_element
        procedure, public :: add_torsional_spring => ldm_add_torsional_spring
        procedure, public :: add_torsional_damper => ldm_add_torsional_damper
        procedure, public :: get_axial_element_count => ldm_get_axial_count
        procedure, public :: get_torsional_element_count => ldm_get_torsional_count
        procedure, public :: get_axial_element_results => ldm_get_axial_results
        procedure, public :: get_torsional_element_results => &
            ldm_get_torsional_results
        procedure, public :: solve => ldm_solve
    end type

    interface linkage_dynamic_model
        module procedure :: ldm_from_serial
        module procedure :: ldm_from_parallel
    end interface

    type linkage_solve_context
        !! Carries model data and optional loads into integrator callbacks.
        class(linkage_dynamic_model), pointer :: model => null()
            !! The dynamic model being integrated.
        real(real64), dimension(3) :: gravity = 0.0d0
            !! The uniform gravitational acceleration in world coordinates.
        real(real64), allocatable, dimension(:,:) :: body_force
            !! Constant world-frame forces, shaped 3-by-number-of-bodies.
        real(real64), allocatable, dimension(:,:) :: body_torque
            !! Constant body-frame torques, shaped 3-by-number-of-bodies.
        integer(int32) :: prescribed_body = 0
            !! The one-based planar body whose absolute angle is prescribed, or
            !! zero when no prescribed motion is active.
        procedure(linkage_prescribed_motion), pointer, nopass :: &
            prescribed_motion => null()
            !! The callback defining the prescribed body's absolute angle, or a
            !! null pointer when no prescribed motion is active.
		class(*), pointer :: user_args
			!! User-specified data to pass along.
    end type

contains
! ------------------------------------------------------------------------------
pure function axial_zero_force(this, length, length_rate) result(rst)
    !! Evaluates the default axial force law, which produces no force.
    class(axial_force_element), intent(in) :: this
        !! The axial element being evaluated.
    real(real64), intent(in) :: length, length_rate
        !! The current element length and its time derivative. They are unused
        !! by the zero-force law.
    real(real64) :: rst
        !! The signed axial force, equal to zero.
    rst = 0.0d0
end function

pure function linear_spring_force(this, length, length_rate) result(rst)
    !! Evaluates the signed force produced by a linear axial spring.
    class(linear_spring), intent(in) :: this
        !! The spring whose stiffness and free length define the force law.
    real(real64), intent(in) :: length, length_rate
        !! The current element length and its time derivative. The rate is
        !! unused by a linear spring.
    real(real64) :: rst
        !! The signed spring force, positive in the element's point-1 to
        !! point-2 direction.
    rst = this%stiffness * (length - this%free_length)
end function

pure function linear_damper_force(this, length, length_rate) result(rst)
    !! Evaluates the signed force produced by a linear axial damper.
    class(linear_damper), intent(in) :: this
        !! The damper whose damping coefficient defines the force law.
    real(real64), intent(in) :: length, length_rate
        !! The current element length and the relative rate of separation of
        !! the two attachment points. The length is unused by this law.
    real(real64) :: rst
        !! The signed damping force, positive in the element's point-1 to
        !! point-2 direction.
    rst = this%damping * length_rate
end function

pure function torsional_zero_torque(this, angle, angle_rate) result(rst)
    !! Evaluates the default torsional law, which produces no torque.
    class(torsional_force_element), intent(in) :: this
        !! The torsional element being evaluated.
    real(real64), intent(in) :: angle, angle_rate
        !! The relative joint angle and angular rate. They are unused by the
        !! zero-torque law.
    real(real64) :: rst
        !! The signed joint torque, equal to zero.
    rst = 0.0d0
end function

pure function torsional_spring_torque(this, angle, angle_rate) result(rst)
    !! Evaluates the restoring torque produced by a linear torsional spring.
    class(torsional_spring), intent(in) :: this
        !! The spring whose stiffness and free angle define the torque law.
    real(real64), intent(in) :: angle, angle_rate
        !! The relative joint angle and angular rate. The rate is unused by a
        !! torsional spring.
    real(real64) :: rst
        !! The signed torque applied to the child body about the joint axis.
    rst = -this%stiffness * (angle - this%free_angle)
end function

pure function torsional_damper_torque(this, angle, angle_rate) result(rst)
    !! Evaluates the resisting torque produced by a linear torsional damper.
    class(torsional_damper), intent(in) :: this
        !! The damper whose damping coefficient defines the torque law.
    real(real64), intent(in) :: angle, angle_rate
        !! The relative joint angle and angular rate. The angle is unused by
        !! this law.
    real(real64) :: rst
        !! The signed torque applied to the child body about the joint axis.
    rst = -this%damping * angle_rate
end function
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
    !! Gets the number of moving rigid bodies in the dynamic model.
    class(linkage_dynamic_model), intent(in) :: this
        !! The dynamic linkage model.
    integer(int32) :: rst
        !! The number of moving bodies; ground is not included.
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
    !! Gets the number of scalar constraints in the dynamic model.
    class(linkage_dynamic_model), intent(in) :: this
        !! The dynamic linkage model.
    integer(int32) :: rst
        !! The total number of joint and planar-body constraints.
    rst = this%m_constraint_count
end function

! ------------------------------------------------------------------------------
function ldm_get_initial_state(this) result(rst)
    !! Gets a copy of the constraint-compatible state used to construct the
    !! dynamic model.
    class(linkage_dynamic_model), intent(in) :: this
        !! The dynamic linkage model.
    type(variational_state) :: rst
        !! A copy of the model's constraint-compatible initial state.
    rst = this%m_initial_state
end function

! ------------------------------------------------------------------------------
subroutine ldm_add_axial_element(this, element)
    !! Adds an extensible axial force element to the dynamic model.
    class(linkage_dynamic_model), intent(inout) :: this
        !! The model that receives the element.
    class(axial_force_element), intent(in) :: element
        !! The axial element to add. Its body indices must identify two
        !! distinct valid bodies or ground.
    type(axial_element_container), allocatable, dimension(:) :: buffer
    integer(int32) :: n

    call validate_body_pair(this, element%body_1, element%body_2)
    n = this%get_axial_element_count()
    allocate(buffer(n + 1))
    if (n > 0) buffer(1:n) = this%m_axial_elements
    allocate(buffer(n + 1)%item, source = element)
    call move_alloc(buffer, this%m_axial_elements)
end subroutine

subroutine ldm_add_linear_spring(this, element)
    !! Adds a validated linear axial spring to the dynamic model.
    class(linkage_dynamic_model), intent(inout) :: this
        !! The model that receives the spring.
    type(linear_spring), intent(in) :: element
        !! The spring to add. Its stiffness and free length must be
        !! nonnegative.
    if (element%stiffness < 0.0d0 .or. element%free_length < 0.0d0) &
        error stop DYN_INVALID_INPUT_ERROR
    call this%add_axial_element(element)
end subroutine

subroutine ldm_add_linear_damper(this, element)
    !! Adds a validated linear axial damper to the dynamic model.
    class(linkage_dynamic_model), intent(inout) :: this
        !! The model that receives the damper.
    type(linear_damper), intent(in) :: element
        !! The damper to add. Its damping coefficient must be nonnegative.
    if (element%damping < 0.0d0) error stop DYN_INVALID_INPUT_ERROR
    call this%add_axial_element(element)
end subroutine

subroutine ldm_add_torsional_element(this, element)
    !! Adds an extensible torsional element bound to a revolute joint.
    class(linkage_dynamic_model), intent(inout) :: this
        !! The model that receives the element.
    class(torsional_force_element), intent(in) :: element
        !! The torsional element to add. Its one-based joint index must refer
        !! to a revolute joint in this model.
    type(torsional_element_container), allocatable, dimension(:) :: buffer
    integer(int32) :: n

    if (element%joint_index < 1 .or. &
        element%joint_index > this%get_joint_count()) &
        error stop DYN_INVALID_INPUT_ERROR
    if (this%m_joints(element%joint_index)%joint_type /= REVOLUTE_JOINT) &
        error stop DYN_INVALID_INPUT_ERROR
    n = this%get_torsional_element_count()
    allocate(buffer(n + 1))
    if (n > 0) buffer(1:n) = this%m_torsional_elements
    allocate(buffer(n + 1)%item, source = element)
    call move_alloc(buffer, this%m_torsional_elements)
end subroutine

subroutine ldm_add_torsional_spring(this, element)
    !! Adds a validated linear torsional spring to the dynamic model.
    class(linkage_dynamic_model), intent(inout) :: this
        !! The model that receives the spring.
    type(torsional_spring), intent(in) :: element
        !! The spring to add. Its stiffness must be nonnegative.
    if (element%stiffness < 0.0d0) error stop DYN_INVALID_INPUT_ERROR
    call this%add_torsional_element(element)
end subroutine

subroutine ldm_add_torsional_damper(this, element)
    !! Adds a validated linear torsional damper to the dynamic model.
    class(linkage_dynamic_model), intent(inout) :: this
        !! The model that receives the damper.
    type(torsional_damper), intent(in) :: element
        !! The damper to add. Its damping coefficient must be nonnegative.
    if (element%damping < 0.0d0) error stop DYN_INVALID_INPUT_ERROR
    call this%add_torsional_element(element)
end subroutine

pure function ldm_get_axial_count(this) result(rst)
    !! Gets the number of axial force elements in the model.
    class(linkage_dynamic_model), intent(in) :: this
        !! The dynamic linkage model.
    integer(int32) :: rst
        !! The number of stored axial elements.
    rst = 0
    if (allocated(this%m_axial_elements)) rst = size(this%m_axial_elements)
end function

pure function ldm_get_torsional_count(this) result(rst)
    !! Gets the number of torsional force elements in the model.
    class(linkage_dynamic_model), intent(in) :: this
        !! The dynamic linkage model.
    integer(int32) :: rst
        !! The number of stored torsional elements.
    rst = 0
    if (allocated(this%m_torsional_elements)) &
        rst = size(this%m_torsional_elements)
end function

function ldm_get_axial_results(this, state) result(rst)
    !! Gets instantaneous lengths, rates, and signed forces for axial elements.
    class(linkage_dynamic_model), intent(in) :: this
        !! The dynamic linkage model containing the elements.
    type(variational_state), intent(in) :: state
        !! The state at which each element is evaluated.
    type(axial_element_result), allocatable, dimension(:) :: rst
        !! One result for each axial element, in insertion order.
    integer(int32) :: i
    real(real64), dimension(3) :: direction, arm_1, arm_2

    allocate(rst(this%get_axial_element_count()))
    do i = 1, size(rst)
        call evaluate_axial_element(this, state, this%m_axial_elements(i)%item, &
            rst(i), direction, arm_1, arm_2)
    end do
end function

function ldm_get_torsional_results(this, state) result(rst)
    !! Gets instantaneous angles, twist rates, and torques for torsional elements.
    class(linkage_dynamic_model), intent(in) :: this
        !! The dynamic linkage model containing the elements.
    type(variational_state), intent(in) :: state
        !! The state at which each element is evaluated.
    type(torsional_element_result), allocatable, dimension(:) :: rst
        !! One result for each torsional element, in insertion order.
    integer(int32) :: i
    real(real64), dimension(3) :: axis

    allocate(rst(this%get_torsional_element_count()))
    do i = 1, size(rst)
        call evaluate_torsional_element(this, state, &
            this%m_torsional_elements(i)%item, rst(i), axis)
    end do
end function

! ------------------------------------------------------------------------------
function ldm_constraint_residual(this, state) result(rst)
    !! Evaluates every joint and planar constraint for a supplied state.
    class(linkage_dynamic_model), intent(in) :: this
        !! The dynamic linkage model whose constraints are evaluated.
    type(variational_state), intent(in) :: state
        !! The body state at which the constraints are evaluated.
    real(real64), allocatable, dimension(:) :: rst
        !! The constraint residual vector in the model's constraint ordering.
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
        !! The planar dynamic linkage model.
    type(variational_state), intent(in) :: state
        !! The state used to orient a prismatic joint's axis.
    type(dynamic_joint), intent(in) :: descriptor
        !! The planar joint descriptor whose reaction is being reconstructed.
    real(real64), intent(in), dimension(:) :: multipliers
        !! The multipliers associated with this joint's planar constraints.
    type(joint_reaction), intent(out) :: reaction
        !! The resulting world-frame force and moment on the child link.

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
        !! The joint kind whose multiplier layout is being decoded.
    real(real64), intent(in), dimension(:) :: multipliers
        !! The spatial joint multipliers in the integrator's ordering.
    real(real64), intent(out), dimension(3) :: force, moment
        !! The force and moment in the parent joint frame.

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
    body_force, body_torque, prescribed_body, prescribed_motion, multipliers, &
	args) result(rst)
    !! Integrates the linkage dynamics under an optional uniform world-frame
    !! gravitational acceleration and optional constant body loads.
    class(linkage_dynamic_model), intent(in), target :: this
        !! The dynamic linkage model to integrate.
    type(variational_integrator), intent(in) :: integrator
        !! The variational integrator used to advance the model.
    real(real64), intent(in) :: dt
        !! The time step in seconds.
    integer(int32), intent(in) :: ntime
        !! The number of time steps to integrate.
    type(variational_state), intent(in), optional :: initial_state
        !! An optional starting state; the model's initial state is used when
        !! this argument is absent.
    real(real64), intent(in), optional, dimension(3) :: gravity
        !! Optional constant gravitational acceleration in world coordinates.
    real(real64), intent(in), optional, dimension(:,:) :: body_force
        !! Constant 3-by-nbody world-frame force array.
    real(real64), intent(in), optional, dimension(:,:) :: body_torque
        !! Constant 3-by-nbody body-frame torque array.
    integer(int32), intent(in), optional :: prescribed_body
        !! The moving body whose absolute planar angle is prescribed.
    procedure(linkage_prescribed_motion), pointer, intent(in), optional :: prescribed_motion
        !! The prescribed absolute planar angle as a function of time.
    real(real64), allocatable, intent(out), optional, dimension(:,:) :: multipliers
        !! Constraint multipliers for each completed time step. When a motion is
        !! prescribed, the last row is the required motor torque.
	class(*), intent(inout), target, optional :: args
		!! A mechanism for the caller to pass information to/from the 
		!! user-defined routines (e.g. presribed_motion).
    type(variational_state), allocatable, dimension(:) :: rst
            !! The state at the initial time followed by the state after each
            !! completed integration step.

    type(linkage_solve_context) :: context
    type(variational_state) :: state
    integer(int32) :: constraint_count
    procedure(variational_force), pointer :: frc
    procedure(variational_constraint), pointer :: constraint

    state = this%m_initial_state
    if (present(initial_state)) state = initial_state
    context%model => this
	context%user_args => null()
	if (present(args)) context%user_args => args
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
    frc => linkage_gravity
    constraint => linkage_constraints
    rst = integrator%solve(this%m_bodies, state, dt, ntime, &
        constraint_count = constraint_count, &
        constraint = constraint, &
        force_function = frc, multipliers = multipliers, &
        args = context)
end function

! ------------------------------------------------------------------------------
subroutine linkage_constraints(state, value, args)
    !! Evaluates the model constraints for the variational integrator.
    type(variational_state), intent(in) :: state
        !! The state at which the constraints are evaluated.
    real(real64), intent(out), dimension(:) :: value
        !! The output constraint residual vector, including an optional
        !! prescribed-motion residual as its final entry.
    class(*), intent(inout), optional :: args
        !! A linkage_solve_context containing the model and optional motion
        !! prescription.

    select type (context => args)
    type is (linkage_solve_context)
        value(1:context%model%m_constraint_count) = &
            context%model%constraint_residual(state)
        if (associated(context%prescribed_motion)) then
			if (associated(context%user_args)) then
				value(size(value)) = planar_body_angle(state, &
					context%prescribed_body) - &
					context%prescribed_motion(state%time, context%user_args)
			else
				value(size(value)) = planar_body_angle(state, &
					context%prescribed_body) - &
					context%prescribed_motion(state%time)
			end if
        end if
    class default
        error stop DYN_INVALID_INPUT_ERROR
    end select
end subroutine

! ------------------------------------------------------------------------------
function planar_body_angle(state, body) result(rst)
    !! Gets the absolute angle of a planar body about the world z-axis.
    type(variational_state), intent(in) :: state
        !! The body state whose orientation is queried.
    integer(int32), intent(in) :: body
        !! The one-based moving-body index.
    real(real64) :: rst
        !! The body's absolute angle in radians, in the range returned by
        !! atan2.
    real(real64), dimension(3,3) :: rotation

    rotation = state%orientation(body)%to_matrix()
    rst = atan2(rotation(2,1), rotation(1,1))
end function

! ------------------------------------------------------------------------------
subroutine linkage_gravity(t, state, force, torque, args)
    !! Builds gravity, applied-load, and force-element loads for one state.
    real(real64), intent(in) :: t
        !! The current simulation time. It is unused by the current load laws.
    type(variational_state), intent(in) :: state
        !! The state used to evaluate force and torque elements.
    real(real64), intent(out), dimension(:,:) :: force, torque
        !! The output world-frame forces and body-frame torques for each body.
    class(*), intent(inout), optional :: args
        !! A linkage_solve_context containing gravity, body loads, and the
        !! dynamic model.
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
        call apply_force_elements(context%model, state, force, torque)
    class default
        error stop DYN_INVALID_INPUT_ERROR
    end select
end subroutine

! ------------------------------------------------------------------------------
subroutine apply_force_elements(model, state, force, torque)
    !! Adds all axial-element forces and torsional-element torques to loads.
    class(linkage_dynamic_model), intent(in) :: model
        !! The dynamic model containing the force elements.
    type(variational_state), intent(in) :: state
        !! The state at which element loads are evaluated.
    real(real64), intent(inout), dimension(:,:) :: force, torque
        !! The accumulated world-frame forces and body-frame torques to update.
    integer(int32) :: i
    real(real64), dimension(3) :: direction, arm_1, arm_2, element_force, axis
    type(axial_element_result) :: axial_result
    type(torsional_element_result) :: torsional_result

    do i = 1, model%get_axial_element_count()
        call evaluate_axial_element(model, state, model%m_axial_elements(i)%item, &
            axial_result, direction, arm_1, arm_2)
        element_force = axial_result%force * direction
        call apply_point_force(model, state, &
            model%m_axial_elements(i)%item%body_1, arm_1, element_force, &
            force, torque)
        call apply_point_force(model, state, &
            model%m_axial_elements(i)%item%body_2, arm_2, -element_force, &
            force, torque)
    end do

    do i = 1, model%get_torsional_element_count()
        call evaluate_torsional_element(model, state, &
            model%m_torsional_elements(i)%item, torsional_result, axis)
        call apply_axis_torque(state, &
            model%m_joints(model%m_torsional_elements(i)%item%joint_index)%child_body, &
            torsional_result%torque * axis, torque)
        call apply_axis_torque(state, &
            model%m_joints(model%m_torsional_elements(i)%item%joint_index)%parent_body, &
            -torsional_result%torque * axis, torque)
    end do
end subroutine

subroutine evaluate_axial_element(model, state, element, result, direction, &
    arm_1, arm_2)
    !! Evaluates an axial element and computes its current geometry.
    class(linkage_dynamic_model), intent(in) :: model
        !! The dynamic model containing the element's bodies.
    type(variational_state), intent(in) :: state
        !! The state at which the attachment points are evaluated.
    class(axial_force_element), intent(in) :: element
        !! The element whose length, rate, and force are computed.
    type(axial_element_result), intent(out) :: result
        !! The current length, length rate, and signed element force.
    real(real64), intent(out), dimension(3) :: direction, arm_1, arm_2
        !! The world direction from point 1 to point 2 and the body-frame
        !! moment arms for the two attachment points.
    real(real64), dimension(3) :: p1, p2, v1, v2

    call attachment_state(model, state, element%body_1, element%point_1, &
        p1, v1, arm_1)
    call attachment_state(model, state, element%body_2, element%point_2, &
        p2, v2, arm_2)
    direction = p2 - p1
    result%length = norm2(direction)
    if (result%length <= sqrt(epsilon(1.0d0))) &
        error stop DYN_INVALID_INPUT_ERROR
    direction = direction / result%length
    result%length_rate = dot_product(v2 - v1, direction)
    result%force = element%evaluate_force(result%length, result%length_rate)
end subroutine

subroutine attachment_state(model, state, body, point, position, velocity, arm)
    !! Transforms an attachment point into world position and velocity.
    class(linkage_dynamic_model), intent(in) :: model
        !! The dynamic model providing body mass-center locations.
    type(variational_state), intent(in) :: state
        !! The state used to transform the point.
    integer(int32), intent(in) :: body
        !! The one-based body index, or zero for a fixed ground point.
    real(real64), intent(in), dimension(3) :: point
        !! The point in body coordinates, or world coordinates for ground.
    real(real64), intent(out), dimension(3) :: position, velocity, arm
        !! The world position and velocity of the point, and its body-frame
        !! moment arm from the body's center of mass.
    real(real64), dimension(3,3) :: rotation
    real(real64), dimension(3) :: omega_world

    if (body == 0) then
        position = point
        velocity = 0.0d0
        arm = 0.0d0
    else
        rotation = state%orientation(body)%to_matrix()
        arm = matmul(rotation, point - model%m_bodies(body)%cg)
        position = state%position(:,body) + arm
        omega_world = matmul(rotation, state%angular_velocity(:,body))
        velocity = state%velocity(:,body) + cross_product(omega_world, arm)
    end if
end subroutine

subroutine apply_point_force(model, state, body, arm, applied, force, torque)
    !! Adds a point force and its moment to one moving body's load arrays.
    class(linkage_dynamic_model), intent(in) :: model
        !! The dynamic model providing body indexing.
    type(variational_state), intent(in) :: state
        !! The state used to rotate the resulting moment into body coordinates.
    integer(int32), intent(in) :: body
        !! The one-based body receiving the load, or zero for ground.
    real(real64), intent(in), dimension(3) :: arm, applied
        !! The body-frame moment arm and world-frame applied force.
    real(real64), intent(inout), dimension(:,:) :: force, torque
        !! The accumulated world-frame force and body-frame torque arrays.
    real(real64), dimension(3,3) :: rotation

    if (body == 0) return
    force(:,body) = force(:,body) + applied
    rotation = state%orientation(body)%to_matrix()
    torque(:,body) = torque(:,body) + &
        matmul(transpose(rotation), cross_product(arm, applied))
end subroutine

subroutine evaluate_torsional_element(model, state, element, result, axis)
    !! Evaluates a torsional element about its associated revolute-joint axis.
    class(linkage_dynamic_model), intent(in) :: model
        !! The dynamic model containing the referenced joint.
    type(variational_state), intent(in) :: state
        !! The state at which the joint angle and rate are evaluated.
    class(torsional_force_element), intent(in) :: element
        !! The element whose torque law is evaluated.
    type(torsional_element_result), intent(out) :: result
        !! The current relative angle, angular rate, and signed torque.
    real(real64), intent(out), dimension(3) :: axis
        !! The joint axis expressed in world coordinates.
    type(dynamic_joint) :: descriptor
    real(real64), dimension(4,4) :: parent, child, relative
    real(real64), dimension(3) :: omega_parent, omega_child

    descriptor = model%m_joints(element%joint_index)
    parent = joint_world_transform(model, state, descriptor%parent_body, &
        descriptor%parent_frame)
    child = joint_world_transform(model, state, descriptor%child_body, &
        descriptor%child_frame)
    relative = matmul(transform_inverse(parent), child)
    result%angle = atan2(relative(2,1), relative(1,1))
    axis = parent(1:3,3)
    call body_angular_velocity(state, descriptor%parent_body, omega_parent)
    call body_angular_velocity(state, descriptor%child_body, omega_child)
    result%angle_rate = dot_product(axis, omega_child - omega_parent)
    result%torque = element%evaluate_torque(result%angle, result%angle_rate)
end subroutine

subroutine body_angular_velocity(state, body, omega)
    !! Gets a moving body's angular velocity in world coordinates.
    type(variational_state), intent(in) :: state
        !! The state containing the body's orientation and angular velocity.
    integer(int32), intent(in) :: body
        !! The one-based body index, or zero for a fixed ground body.
    real(real64), intent(out), dimension(3) :: omega
        !! The body's world-frame angular velocity.
    if (body == 0) then
        omega = 0.0d0
    else
        omega = matmul(state%orientation(body)%to_matrix(), &
            state%angular_velocity(:,body))
    end if
end subroutine

subroutine apply_axis_torque(state, body, applied, torque)
    !! Adds a world-frame axis torque to one moving body's torque array.
    type(variational_state), intent(in) :: state
        !! The state used to rotate the torque into body coordinates.
    integer(int32), intent(in) :: body
        !! The one-based body receiving the torque, or zero for ground.
    real(real64), intent(in), dimension(3) :: applied
        !! The applied torque in world coordinates.
    real(real64), intent(inout), dimension(:,:) :: torque
        !! The accumulated body-frame torque array to update.
    if (body == 0) return
    torque(:,body) = torque(:,body) + &
        matmul(transpose(state%orientation(body)%to_matrix()), applied)
end subroutine

subroutine validate_body_pair(model, body_1, body_2)
    !! Validates two distinct body or ground indices for an axial element.
    class(linkage_dynamic_model), intent(in) :: model
        !! The model defining the valid moving-body index range.
    integer(int32), intent(in) :: body_1, body_2
        !! The two indices to validate; zero denotes ground.
    if (body_1 < 0 .or. body_1 > model%get_body_count() .or. &
        body_2 < 0 .or. body_2 > model%get_body_count() .or. &
        body_1 == body_2) error stop DYN_INVALID_INPUT_ERROR
end subroutine

! ------------------------------------------------------------------------------
subroutine append_planar_body_constraints(model, state, body, value, index)
    !! Appends position and orientation constraints for one planar body.
    class(linkage_dynamic_model), intent(in) :: model
        !! The planar dynamic model whose body pose is evaluated.
    type(variational_state), intent(in) :: state
        !! The state at which the body constraints are evaluated.
    integer(int32), intent(in) :: body
        !! The one-based moving-body index.
    real(real64), intent(inout), dimension(:) :: value
        !! The constraint vector receiving the three scalar residuals.
    integer(int32), intent(inout) :: index
        !! The one-based insertion position, advanced past the appended values.
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
    !! Appends the residuals for one planar or spatial joint.
    class(linkage_dynamic_model), intent(in) :: model
        !! The dynamic model whose planar/spatial mode is selected.
    type(variational_state), intent(in) :: state
        !! The state at which the joint residual is evaluated.
    type(dynamic_joint), intent(in) :: descriptor
        !! The joint descriptor defining type, bodies, and frames.
    real(real64), intent(inout), dimension(:) :: value
        !! The constraint vector receiving the joint residuals.
    integer(int32), intent(inout) :: index
        !! The one-based insertion position, advanced past the appended values.
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
    !! Appends spatial joint residuals from a relative joint transform.
    integer(int32), intent(in) :: joint_type
        !! The joint kind determining the residual components and count.
    real(real64), intent(in), dimension(4,4) :: relative
        !! The child joint frame relative to the parent joint frame.
    real(real64), intent(inout), dimension(:) :: value
        !! The constraint vector receiving the spatial residuals.
    integer(int32), intent(inout) :: index
        !! The one-based insertion position, advanced past the appended values.
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
    !! Converts a body-local or ground joint frame into world coordinates.
    class(linkage_dynamic_model), intent(in) :: model
        !! The dynamic model providing body center-of-mass offsets.
    type(variational_state), intent(in) :: state
        !! The state providing the selected body's pose.
    integer(int32), intent(in) :: body
        !! The one-based moving-body index, or zero for ground.
    real(real64), intent(in), dimension(4,4) :: frame
        !! The joint frame relative to the selected body, or already in world
        !! coordinates when body is zero.
    real(real64), dimension(4,4) :: rst, transform
        !! The resulting world-frame homogeneous transform.
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
    !! Gets a moving body's orientation and reference-frame origin in world
    !! coordinates.
    class(linkage_dynamic_model), intent(in) :: model
        !! The dynamic model providing the body's center-of-mass location.
    type(variational_state), intent(in) :: state
        !! The state containing the body's position and orientation.
    integer(int32), intent(in) :: body
        !! The one-based moving-body index.
    real(real64), intent(out), dimension(3,3) :: rotation
        !! The body-to-world rotation matrix.
    real(real64), intent(out), dimension(3) :: origin
        !! The world position of the body's local coordinate origin.

    rotation = state%orientation(body)%to_matrix()
    origin = state%position(:,body) - &
        matmul(rotation, model%m_bodies(body)%cg)
end subroutine

! ------------------------------------------------------------------------------
subroutine allocate_state(state, nbody)
    !! Allocates and initializes a variational state for moving bodies.
    type(variational_state), intent(out) :: state
        !! The state to allocate and initialize.
    integer(int32), intent(in) :: nbody
        !! The number of moving bodies represented by the state.
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
    !! Initializes one body's position and orientation from a homogeneous pose.
    type(variational_state), intent(inout) :: state
        !! The state receiving the body's pose.
    integer(int32), intent(in) :: body
        !! The one-based body index to initialize.
    class(link), intent(in) :: properties
        !! The link mass properties, including its center-of-mass offset.
    real(real64), intent(in), dimension(4,4) :: transform
        !! The world pose of the body's local coordinate origin.

    state%orientation(body) = quaternion(transform(1:3,1:3))
    state%position(:,body) = transform(1:3,4) + &
        matmul(transform(1:3,1:3), properties%cg)
end subroutine

! ------------------------------------------------------------------------------
pure subroutine copy_mass_properties(body, source)
    !! Copies the inherited rigid-body properties from a linkage link.
    type(rigid_body), intent(out) :: body
        !! The rigid body receiving mass, center-of-mass, and inertia data.
    class(link), intent(in) :: source
        !! The linkage link supplying those properties.

    body%mass = source%mass
    body%cg = source%cg
    body%inertia = source%inertia
end subroutine

! ------------------------------------------------------------------------------
pure function serial_joint_motion(joint_type, q) result(rst)
    !! Builds the homogeneous motion transform for one serial joint variable.
    integer(int32), intent(in) :: joint_type
        !! The joint kind; revolute and prismatic joints are supported.
    real(real64), intent(in) :: q
        !! The joint variable, in radians for a revolute joint or length units
        !! for a prismatic joint.
    real(real64), dimension(4,4) :: rst
        !! The joint motion transform.
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
    !! Sums the scalar constraint counts for a collection of joints.
    type(dynamic_joint), intent(in), dimension(:) :: joints
        !! The joint descriptors whose constraints are counted.
    logical, intent(in) :: planar
        !! Whether planar rather than spatial constraint counts are required.
    integer(int32) :: rst
        !! The total number of scalar joint constraints.
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
        !! The joint kind whose scalar constraint count is requested.
    logical, intent(in) :: planar
        !! Whether to use planar or spatial joint constraints.
    integer(int32) :: rst
        !! The number of scalar constraints, or an error for an unsupported
        !! joint/mode combination.

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

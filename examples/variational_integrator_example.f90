! This example simulates a planar double pendulum using maximal coordinates.
! Both connecting rods are rigid bodies with distributed mass and rotational
! inertia. Six holonomic constraints pin the first rod to ground and connect
! the distal end of the first rod to the proximal end of the second rod.

program example
	use iso_fortran_env, only : int32, real64
	use dynamics
	use fplot_core
	implicit none

	type double_pendulum_parameters
		!! Defines the geometry and gravity used by the example callbacks.
		real(real64) :: length1
			!! The length of the first rod.
		real(real64) :: length2
			!! The length of the second rod.
		real(real64) :: gravity
			!! The gravitational acceleration magnitude.
	end type

	! Model Parameters
	integer(int32), parameter :: nbody = 2
	integer(int32), parameter :: nconstraint = 6
	integer(int32), parameter :: ntime = 801
	real(real64), parameter :: dt = 5.0d-3
	real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
	real(real64), parameter :: length1 = 1.0d0
	real(real64), parameter :: length2 = 0.8d0
	real(real64), parameter :: width1 = 0.06d0
	real(real64), parameter :: width2 = 0.05d0
	real(real64), parameter :: mass1 = 1.5d0
	real(real64), parameter :: mass2 = 1.0d0
	real(real64), parameter :: gravity = 9.80665d0
	real(real64), parameter :: angle1_initial = 35.0d0 * pi / 180.0d0
	real(real64), parameter :: angle2_initial = -20.0d0 * pi / 180.0d0

	! Local Variables
	integer(int32) :: i
	real(real64), dimension(3,3) :: inertia1, inertia2, rotation1, rotation2
	real(real64), dimension(3) :: direction1, direction2
	real(real64), allocatable, dimension(:) :: time, angle1, angle2
	type(rigid_body), dimension(nbody) :: bodies
	type(variational_state) :: initial_state
	type(variational_state), allocatable, dimension(:) :: solution
	type(variational_integrator) :: integrator
	type(double_pendulum_parameters) :: parameters
	procedure(variational_constraint), pointer :: constraint_ptr
	procedure(variational_force), pointer :: force_ptr
	procedure(variational_constraint_jacobian), pointer :: jacobian_ptr
	type(plot_2d) :: plt
    class(plot_axis), pointer :: xAxis
	type(legend), pointer :: lgnd

	! Define rectangular-prism inertia tensors about each rod's center of mass.
	! The body y-axis lies along the rod, and the x- and z-dimensions are equal
	! to the specified width.
	inertia1 = rod_inertia(mass1, length1, width1)
	inertia2 = rod_inertia(mass2, length2, width2)
	bodies(1) = rigid_body(mass1, inertia1)
	bodies(2) = rigid_body(mass2, inertia2)

	! Establish a compatible initial configuration. At zero angle each rod
	! points vertically downward along its negative body y-axis.
	call initialize_variational_state(initial_state, nbody)
	initial_state%orientation(1) = quaternion(angle1_initial, &
		[0.0d0, 0.0d0, 1.0d0])
	initial_state%orientation(2) = quaternion(angle2_initial, &
		[0.0d0, 0.0d0, 1.0d0])
	direction1 = [sin(angle1_initial), -cos(angle1_initial), 0.0d0]
	direction2 = [sin(angle2_initial), -cos(angle2_initial), 0.0d0]
	initial_state%position(:,1) = 0.5d0 * length1 * direction1
	initial_state%position(:,2) = length1 * direction1 + &
		0.5d0 * length2 * direction2

	! Select the graph-factorized Newton solver from the paper. The force and
	! constraint callbacks receive the model parameters through args.
	parameters = double_pendulum_parameters(length1, length2, gravity)
	integrator%settings%linear_solver = VI_GRAPH_FACTORIZED_SOLVER
	constraint_ptr => pendulum_constraints
	force_ptr => gravity_forces
	jacobian_ptr => pendulum_constraint_jacobian
	solution = integrator%solve(bodies, initial_state, dt, ntime, &
		constraint_count = nconstraint, &
		constraint = constraint_ptr, &
		force_function = force_ptr, &
		constraint_jacobian = jacobian_ptr, &
		args = parameters)

	! Recover the two planar angles from the body-to-world rotation matrices.
	allocate(time(ntime), angle1(ntime), angle2(ntime))
	do i = 1, ntime
		time(i) = solution(i)%time
		rotation1 = solution(i)%orientation(1)%to_matrix()
		rotation2 = solution(i)%orientation(2)%to_matrix()
		angle1(i) = atan2(-rotation1(1,2), rotation1(2,2)) * 1.8d2 / pi
		angle2(i) = atan2(-rotation2(1,2), rotation2(2,2)) * 1.8d2 / pi
	end do

	print "(A,F8.3,A)", "Simulated ", time(ntime), " seconds."
	print "(A,2F10.3)", "Final rod angles [deg]: ", &
		angle1(ntime), angle2(ntime)

	! Plot the absolute angle of each massive rod as a function of time.
	call plt%initialize()
	call plt%set_x_axis_title("Time [s]")
	call plt%set_y_axis_title("Rod Angle [deg]")
	lgnd => plt%get_legend()
	call lgnd%set_is_visible(.true.)
	call lgnd%set_draw_border(.false.)
    xAxis => plt%get_x_axis()
    call xAxis%set_zero_axis(.true.)
	call plt%push(time, angle1, name = "Rod 1")
	call plt%push(time, angle2, name = "Rod 2")
	call plt%draw()

contains
! ------------------------------------------------------------------------------
	pure function rod_inertia(mass, length, width) result(rst)
		!! Computes the center-of-mass inertia tensor of a square-section rod
		!! whose longitudinal axis is the body y-axis.
		real(real64), intent(in) :: mass
			!! The rod mass.
		real(real64), intent(in) :: length
			!! The rod length.
		real(real64), intent(in) :: width
			!! The rod width and depth.
		real(real64), dimension(3,3) :: rst
			!! The body-frame inertia tensor.

		rst = 0.0d0
		rst(1,1) = mass * (length**2 + width**2) / 12.0d0
		rst(2,2) = mass * width**2 / 6.0d0
		rst(3,3) = rst(1,1)
	end function

! ------------------------------------------------------------------------------
	subroutine gravity_forces(t, state, force, torque, args)
		!! Applies gravity at the center of mass of each rod. Consequently,
		!! gravity contributes no direct body-frame torque.
		real(real64), intent(in) :: t
			!! The current time; unused because gravity is constant.
		type(variational_state), intent(in) :: state
			!! The current state; unused because gravity is uniform.
		real(real64), intent(out), dimension(:,:) :: force
			!! The world-frame forces acting on the bodies.
		real(real64), intent(out), dimension(:,:) :: torque
			!! The body-frame torques acting on the bodies.
		class(*), intent(inout), optional :: args
			!! The double-pendulum parameters.

		real(real64) :: acceleration

		select type (model => args)
		type is (double_pendulum_parameters)
			acceleration = model%gravity
		class default
			error stop "Invalid double-pendulum callback data."
		end select

		force = 0.0d0
		torque = 0.0d0
		force(2,1) = -bodies(1)%mass * acceleration
		force(2,2) = -bodies(2)%mass * acceleration
	end subroutine

! ------------------------------------------------------------------------------
	subroutine pendulum_constraints(state, value, args)
		!! Pins the proximal endpoint of rod 1 to the origin and joins the
		!! distal endpoint of rod 1 to the proximal endpoint of rod 2.
		type(variational_state), intent(in) :: state
			!! The maximal-coordinate state to evaluate.
		real(real64), intent(out), dimension(:) :: value
			!! The six endpoint-compatibility residuals.
		class(*), intent(inout), optional :: args
			!! The double-pendulum parameters.

		real(real64) :: first_length, second_length
		real(real64), dimension(3,3) :: first_rotation, second_rotation
		real(real64), dimension(3) :: first_proximal, first_distal, &
			second_proximal

		select type (model => args)
		type is (double_pendulum_parameters)
			first_length = model%length1
			second_length = model%length2
		class default
			error stop "Invalid double-pendulum callback data."
		end select

		first_rotation = state%orientation(1)%to_matrix()
		second_rotation = state%orientation(2)%to_matrix()
		first_proximal = state%position(:,1) + matmul(first_rotation, &
			[0.0d0, 0.5d0 * first_length, 0.0d0])
		first_distal = state%position(:,1) + matmul(first_rotation, &
			[0.0d0, -0.5d0 * first_length, 0.0d0])
		second_proximal = state%position(:,2) + matmul(second_rotation, &
			[0.0d0, 0.5d0 * second_length, 0.0d0])

		value(1:3) = first_proximal
		value(4:6) = first_distal - second_proximal
	end subroutine

! ------------------------------------------------------------------------------
	subroutine pendulum_constraint_jacobian(state, jacobian, args)
		!! Computes the reduced endpoint-constraint Jacobian analytically. For
		!! a body point r and local quaternion-vector variation epsilon,
		!!
		!! $$\frac{\partial R(q)r}{\partial\epsilon}=-2R(q)[r]_{\times}.$$
		type(variational_state), intent(in) :: state
			!! The maximal-coordinate state at which to evaluate the Jacobian.
		real(real64), intent(out), dimension(:,:) :: jacobian
			!! The 6-by-12 reduced constraint Jacobian.
		class(*), intent(inout), optional :: args
			!! The double-pendulum parameters.

		real(real64) :: first_length, second_length
		real(real64), dimension(3,3) :: identity, first_rotation, &
			second_rotation
		real(real64), dimension(3) :: first_proximal, first_distal, &
			second_proximal

		select type (model => args)
		type is (double_pendulum_parameters)
			first_length = model%length1
			second_length = model%length2
		class default
			error stop "Invalid double-pendulum callback data."
		end select

		identity = 0.0d0
		identity(1,1) = 1.0d0
		identity(2,2) = 1.0d0
		identity(3,3) = 1.0d0
		first_rotation = state%orientation(1)%to_matrix()
		second_rotation = state%orientation(2)%to_matrix()
		first_proximal = [0.0d0, 0.5d0 * first_length, 0.0d0]
		first_distal = [0.0d0, -0.5d0 * first_length, 0.0d0]
		second_proximal = [0.0d0, 0.5d0 * second_length, 0.0d0]

		jacobian = 0.0d0
		jacobian(1:3,1:3) = identity
		jacobian(1:3,4:6) = -2.0d0 * matmul(first_rotation, &
			to_skew_symmetric(first_proximal))
		jacobian(4:6,1:3) = identity
		jacobian(4:6,4:6) = -2.0d0 * matmul(first_rotation, &
			to_skew_symmetric(first_distal))
		jacobian(4:6,7:9) = -identity
		jacobian(4:6,10:12) = 2.0d0 * matmul(second_rotation, &
			to_skew_symmetric(second_proximal))
	end subroutine

end program

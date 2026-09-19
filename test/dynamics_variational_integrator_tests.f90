module dynamics_variational_integrator_tests
    use iso_fortran_env, only : real64
    use fortran_test_helper
    use dynamics
    implicit none

    real(real64), save, dimension(3) :: scaled_target_position

contains
! ------------------------------------------------------------------------------
function test_variational_free_body() result(rst)
    !! Tests momentum preservation and quaternion advancement for a free body.
    logical :: rst
        !! True if the free-body update is correct; else, false.
    type(rigid_body) :: bodies(1)
        !! The body properties.
    type(variational_state) :: state
        !! The maximal-coordinate body state.
    type(variational_integrator) :: integrator
        !! The integrator under test.
    real(real64), parameter :: dt = 0.01d0
        !! The integration time step.
    real(real64) :: expected_q(4)
        !! The expected orientation quaternion.

    ! Initialization
    rst = .true.
    bodies(1) = rigid_body(2.0d0)
    call initialize_variational_state(state, 1)
    state%velocity(:,1) = [1.0d0, -2.0d0, 0.5d0]
    state%angular_velocity(:,1) = [0.0d0, 0.0d0, 1.0d0]

    ! Advance one free-body step.
    call integrator%step(bodies, state, dt)
    expected_q = [sqrt(1.0d0 - (0.5d0*dt)**2), 0.0d0, 0.0d0, 0.5d0*dt]

    ! Test the translational and rotational updates.
    if (.not.assert(state%position(:,1), &
        [0.01d0, -0.02d0, 0.005d0], 1.0d-12)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_variational_free_body - position"
    end if
    if (.not.assert(state%velocity(:,1), &
        [1.0d0, -2.0d0, 0.5d0], 1.0d-12)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_variational_free_body - velocity"
    end if
    if (.not.assert(state%orientation(1)%to_array(), expected_q, &
        1.0d-12)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_variational_free_body - orientation"
    end if
end function

! ------------------------------------------------------------------------------
function test_variational_applied_force() result(rst)
    !! Tests the translational discrete momentum balance under constant force.
    logical :: rst
        !! True if the forced update is correct; else, false.
    type(rigid_body) :: bodies(1)
        !! The body properties.
    type(variational_state) :: state
        !! The maximal-coordinate body state.
    type(variational_integrator) :: integrator
        !! The integrator under test.
    procedure(variational_force), pointer :: force_ptr
        !! The applied-force callback under test.

    ! Initialization
    rst = .true.
    bodies(1) = rigid_body(2.0d0)
    call initialize_variational_state(state, 1)
    force_ptr => constant_force

    ! Advance one step under a constant downward force.
    call integrator%step(bodies, state, 0.01d0, &
        force_function = force_ptr)

    ! Test the resulting velocity and position.
    if (.not.assert(state%velocity(:,1), &
        [0.0d0, 0.0d0, -0.0981d0], 1.0d-10)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_variational_applied_force - velocity"
    end if
    if (.not.assert(state%position(:,1), &
        [0.0d0, 0.0d0, -0.000981d0], 1.0d-12)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_variational_applied_force - position"
    end if
end function

! ------------------------------------------------------------------------------
function test_variational_position_constraint() result(rst)
    !! Tests position-level equality-constraint enforcement and multipliers.
    logical :: rst
        !! True if the constrained update is correct; else, false.
    type(rigid_body) :: bodies(1)
        !! The body properties.
    type(variational_state) :: state
        !! The maximal-coordinate body state.
    type(variational_integrator) :: integrator
        !! The integrator under test.
    real(real64), allocatable :: multipliers(:)
        !! The computed equality-constraint multipliers.
    procedure(variational_constraint), pointer :: constraint_ptr
        !! The center-of-mass position constraint under test.

    ! Initialization
    rst = .true.
    bodies(1) = rigid_body(2.0d0)
    call initialize_variational_state(state, 1)
    state%velocity(:,1) = [1.0d0, -2.0d0, 0.5d0]
    constraint_ptr => fixed_position

    ! Advance while fixing all three center-of-mass coordinates.
    call integrator%step(bodies, state, 0.01d0, 3, &
        constraint_ptr, multipliers = multipliers)

    ! Test the constrained state and reaction multipliers.
    if (.not.assert(state%position(:,1), [0.0d0, 0.0d0, 0.0d0], &
        1.0d-10)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_variational_position_constraint - position"
    end if
    if (.not.assert(state%velocity(:,1), [0.0d0, 0.0d0, 0.0d0], &
        1.0d-10)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_variational_position_constraint - velocity"
    end if
    if (.not.assert(multipliers, [-200.0d0, 400.0d0, -100.0d0], &
        1.0d-4)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_variational_position_constraint - multiplier"
    end if
end function

! ------------------------------------------------------------------------------
function test_variational_analytic_constraint_jacobian() result(rst)
    !! Verifies the optional analytic constraint Jacobian path.
    logical :: rst
    type(rigid_body) :: bodies(1)
    type(variational_state) :: state
    type(variational_integrator) :: integrator
    real(real64), allocatable, dimension(:) :: multipliers
    procedure(variational_constraint), pointer :: constraint_ptr
    procedure(variational_constraint_jacobian), pointer :: jacobian_ptr

    rst = .true.
    bodies(1) = rigid_body(2.0d0)
    call initialize_variational_state(state, 1)
    state%velocity(:,1) = [1.0d0, -2.0d0, 0.5d0]
    constraint_ptr => fixed_position
    jacobian_ptr => fixed_position_jacobian
    call integrator%step(bodies, state, 0.01d0, 3, constraint_ptr, &
        constraint_jacobian = jacobian_ptr, &
        multipliers = multipliers)
    if (.not.assert(state%position(:,1), [0.0d0, 0.0d0, 0.0d0], &
        1.0d-10)) then
        rst = .false.
        print "(A)", &
            "TEST FAILED: test_variational_analytic_constraint_jacobian"
    end if
end function

! ------------------------------------------------------------------------------
function test_variational_graph_solver() result(rst)
    !! Tests the graph-factorized Newton solver against the dense reference
    !! solver for a two-body, three-equation relative-position constraint.
    logical :: rst
        !! True if the graph and dense solutions agree; else, false.
    type(rigid_body) :: bodies(2)
        !! The two rigid bodies in the test problem.
    type(variational_state) :: dense_state, graph_state
        !! Identical initial states advanced by the two solver variants.
    type(variational_integrator) :: dense_integrator, graph_integrator
        !! The dense and graph-factorized integrators.
    real(real64), allocatable :: dense_multipliers(:), graph_multipliers(:)
        !! Constraint multipliers returned by each solver.
    procedure(variational_constraint), pointer :: constraint_ptr
        !! The relative-position constraint under test.

    ! Initialization
    rst = .true.
    bodies(1) = rigid_body(1.0d0)
    bodies(2) = rigid_body(1.0d0)
    call initialize_variational_state(dense_state, 2)
    dense_state%position(:,2) = [1.0d0, 0.0d0, 0.0d0]
    dense_state%velocity(:,1) = [1.0d0, 0.0d0, 0.0d0]
    dense_state%velocity(:,2) = [-1.0d0, 0.0d0, 0.0d0]
    graph_state = dense_state
    graph_integrator%settings%linear_solver = VI_GRAPH_FACTORIZED_SOLVER
    constraint_ptr => relative_position_constraint

    ! Advance both copies of the problem.
    call dense_integrator%step(bodies, dense_state, 0.01d0, 3, &
        constraint_ptr, multipliers = dense_multipliers)
    call graph_integrator%step(bodies, graph_state, 0.01d0, 3, &
        constraint_ptr, multipliers = graph_multipliers)

    ! The factorization strategy must not alter the nonlinear solution.
    if (.not.assert(graph_state%position, dense_state%position, 1.0d-10)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_variational_graph_solver - position"
    end if
    if (.not.assert(graph_state%velocity, dense_state%velocity, 1.0d-10)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_variational_graph_solver - velocity"
    end if
    if (.not.assert(graph_multipliers, dense_multipliers, 1.0d-7)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_variational_graph_solver - multiplier"
    end if
end function

! ------------------------------------------------------------------------------
function test_variational_multiplier_history() result(rst)
    !! Verifies that solve returns one multiplier vector for every state,
    !! including the initial point and a look-ahead value at the final point.
    logical :: rst
    type(rigid_body) :: bodies(1)
    type(variational_state) :: state
    type(variational_state), allocatable, dimension(:) :: solution
    type(variational_integrator) :: integrator
    real(real64), allocatable, dimension(:,:) :: multipliers
    procedure(variational_constraint), pointer :: constraint_ptr
    procedure(variational_force), pointer :: force_ptr

    rst = .true.
    bodies(1) = rigid_body(2.0d0)
    call initialize_variational_state(state, 1)
    constraint_ptr => fixed_position
    force_ptr => constant_force
    solution = integrator%solve(bodies, state, 0.01d0, 3, &
        constraint_count = 3, constraint = constraint_ptr, &
        force_function = force_ptr, multipliers = multipliers)
    if (.not.assert(size(multipliers,1), 3)) rst = .false.
    if (.not.assert(size(multipliers,2), size(solution))) rst = .false.
    if (.not.assert(multipliers(:,1), [0.0d0, 0.0d0, 19.62d0], &
        1.0d-6)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_variational_multiplier_history - initial"
    end if
    if (.not.assert(multipliers(:,3), multipliers(:,2), 1.0d-6)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_variational_multiplier_history - final"
    end if
end function

! ------------------------------------------------------------------------------
function test_scaled_constraint_differences() result(rst)
    !! Verifies that translation perturbations remain resolvable at large world
    !! coordinates while rotation perturbations use their independent scale.
    logical :: rst
    type(rigid_body) :: bodies(1)
    type(variational_state) :: state
    type(variational_integrator) :: integrator
    real(real64), allocatable, dimension(:) :: multipliers
    procedure(variational_constraint), pointer :: constraint_ptr

    rst = .true.
    bodies(1) = rigid_body(1.0d0)
    call initialize_variational_state(state, 1)
    scaled_target_position = [1.0d9, -1.0d9, 5.0d8]
    state%position(:,1) = scaled_target_position
    state%velocity(:,1) = [1.0d0, -2.0d0, 0.5d0]
    state%angular_velocity(:,1) = [0.1d0, -0.2d0, 0.3d0]
    integrator%settings%constraint_translation_scale = 1.0d3
    integrator%settings%constraint_rotation_scale = 0.25d0
    constraint_ptr => scaled_fixed_pose
    call integrator%step(bodies, state, 1.0d-3, 6, &
        constraint_ptr, multipliers = multipliers)
    if (maxval(abs(state%position(:,1) - scaled_target_position)) > &
        1.0d-6 .or. norm2(aimag(state%orientation(1))) > 1.0d-8) then
        rst = .false.
        print "(A)", "TEST FAILED: test_scaled_constraint_differences"
    end if
end function

! ------------------------------------------------------------------------------
subroutine constant_force(t, state, force, torque, args)
    !! Supplies the constant force used by the applied-force test.
    real(real64), intent(in) :: t
        !! The current time; unused by this callback.
    type(variational_state), intent(in) :: state
        !! The current state; unused by this callback.
    real(real64), intent(out) :: force(:,:), torque(:,:)
        !! The output force and torque arrays.
    class(*), intent(inout), optional :: args
        !! Optional user data; unused by this callback.

    force = 0.0d0
    torque = 0.0d0
    force(3,1) = -19.62d0
end subroutine

! ------------------------------------------------------------------------------
subroutine fixed_position(state, value, args)
    !! Fixes the center-of-mass position of the first body at the origin.
    type(variational_state), intent(in) :: state
        !! The state at which to evaluate the constraint.
    real(real64), intent(out) :: value(:)
        !! The three position residuals.
    class(*), intent(inout), optional :: args
        !! Optional user data; unused by this callback.

    value = state%position(:,1)
end subroutine

! ------------------------------------------------------------------------------
subroutine fixed_position_jacobian(state, jacobian, args)
    !! Supplies the exact reduced Jacobian for fixed center-of-mass position.
    type(variational_state), intent(in) :: state
    real(real64), intent(out), dimension(:,:) :: jacobian
    class(*), intent(inout), optional :: args

    jacobian = 0.0d0
    jacobian(1,1) = 1.0d0
    jacobian(2,2) = 1.0d0
    jacobian(3,3) = 1.0d0
end subroutine

! ------------------------------------------------------------------------------
subroutine relative_position_constraint(state, value, args)
    !! Evaluates a fixed relative-position constraint between two bodies.
    type(variational_state), intent(in) :: state
        !! The maximal-coordinate state to evaluate.
    real(real64), intent(out) :: value(:)
        !! The three relative-position constraint residuals.
    class(*), intent(inout), optional :: args
        !! Optional user data; unused by this test callback.

    value = state%position(:,2) - state%position(:,1) - &
        [1.0d0, 0.0d0, 0.0d0]
end subroutine

! ------------------------------------------------------------------------------
subroutine scaled_fixed_pose(state, value, args)
    !! Fixes a body pose at large translational coordinates.
    type(variational_state), intent(in) :: state
    real(real64), intent(out), dimension(:) :: value
    class(*), intent(inout), optional :: args

    value(1:3) = state%position(:,1) - scaled_target_position
    value(4:6) = aimag(state%orientation(1))
end subroutine

end module
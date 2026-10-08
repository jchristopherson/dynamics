! Copyright (c) 2022-2026 Jason Christopherson
! SPDX-License-Identifier: MIT
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The Software is provided "as is", without warranty of any kind, express or
! implied, including but not limited to the warranties of merchantability,
! fitness for a particular purpose and noninfringement.
module dynamics_variational_integrator_tests
    use iso_fortran_env, only : int32, real64
    use fortran_test_helper
    use dynamics
    implicit none

    real(real64), save, dimension(3) :: scaled_target_position

contains
! ------------------------------------------------------------------------------
function test_variational_startup() result(rst)
    logical :: rst
    type(rigid_body) :: bodies(1)
    type(variational_state) :: state, continued_state
    type(variational_state), allocatable :: solution(:), continuation(:)
    type(variational_integrator) :: integrator
    type(variational_integrator_info) :: info
    procedure(variational_constraint), pointer :: constraint_ptr
    procedure(variational_force), pointer :: force_ptr
    real(real64), allocatable :: multipliers(:,:)
    real(real64), parameter :: dt = 0.01d0
    real(real64) :: inertia(3,3), physical_momentum(3), discrete_momentum(3), omega(3)
    integer(int32) :: solver, mode

    rst = .true.
    inertia = 0.0d0
    inertia(1,1) = 1.0d0
    inertia(2,2) = 1.0d0
    inertia(3,3) = 1.0d0
    bodies(1) = rigid_body(2.0d0, inertia = inertia)
    constraint_ptr => accelerating_pose
    do solver = VI_DENSE_SOLVER, VI_GRAPH_FACTORIZED_SOLVER
        integrator%settings%linear_solver = solver
        do mode = VI_FORCE_LEFT_ENDPOINT, VI_FORCE_MIDPOINT
            integrator%settings%force_evaluation = mode
            call initialize_variational_state(state, 1)
            continued_state = state
            solution = integrator%solve(bodies, state, dt, 3, &
                constraint_count = 2, constraint = constraint_ptr, &
                multipliers = multipliers, info = info)
            if (.not.info%converged) then
                rst = .false.
                print "(A)", "TEST FAILED: test_variational_startup - convergence"
                return
            end if
            if (.not.assert(multipliers(1,:), [4.0d0, 4.0d0, 4.0d0], 1.0d-6)) rst = .false.
            if (.not.assert(multipliers(2,:), [2.0d0, 2.0d0, 2.0d0], 1.0d-5)) rst = .false.
            if (.not.assert(solution(2)%velocity(1,1), dt, 1.0d-8)) rst = .false.
            if (solution(1)%discrete_momenta .or. .not.state%discrete_momenta) rst = .false.
            continuation = integrator%solve(bodies, continued_state, dt, 2, &
                constraint_count = 2, constraint = constraint_ptr)
            continuation = integrator%solve(bodies, continued_state, dt, 2, &
                constraint_count = 2, constraint = constraint_ptr)
            if (.not.assert(continued_state%position, state%position, 1.0d-8)) rst = .false.
            if (.not.assert(continued_state%angular_velocity, state%angular_velocity, 1.0d-8)) rst = .false.
        end do
    end do

    force_ptr => constant_force
    call initialize_variational_state(state, 1)
    solution = integrator%solve(bodies, state, dt, 3, force_function = force_ptr)
    if (.not.assert(state%position(3,1), -0.5d0 * 9.81d0 * (2.0d0 * dt)**2, 1.0d-9)) rst = .false.
    if (.not.assert(solution(2)%velocity(3,1), -0.5d0 * 9.81d0 * dt, 1.0d-9)) rst = .false.

    call initialize_variational_state(state, 1)
    call integrator%step(bodies, state, dt, force_function = force_ptr, initialize_momenta = .true.)
    if (.not.assert(state%velocity(3,1), -0.5d0 * 9.81d0 * dt, 1.0d-9)) rst = .false.

    inertia(2,2) = 2.0d0
    inertia(3,3) = 3.0d0
    bodies(1) = rigid_body(2.0d0, inertia = inertia)
    call initialize_variational_state(state, 1)
    state%angular_velocity(:,1) = [0.3d0, -0.4d0, 0.5d0]
    physical_momentum = matmul(inertia, state%angular_velocity(:,1))
    call integrator%step(bodies, state, dt, initialize_momenta = .true.)
    omega = state%angular_velocity(:,1)
    discrete_momentum = 0.5d0 * dt * (&
        sqrt(4.0d0 / dt**2 - dot_product(omega, omega)) * matmul(inertia, omega) + &
        cross_product(omega, matmul(inertia, omega)))
    if (.not.assert(discrete_momentum, physical_momentum, 1.0d-9)) rst = .false.

    call initialize_variational_state(state, 1)
    integrator%settings%maximum_iterations = 1
    integrator%settings%tolerance = 1.0d-30
    call integrator%step(bodies, state, dt, force_function = force_ptr, &
        initialize_momenta = .true., info = info)
    if (info%converged .or. state%discrete_momenta .or. state%time /= 0.0d0) rst = .false.
    if (.not.rst) print "(A)", "TEST FAILED: test_variational_startup"
end function

! ------------------------------------------------------------------------------
subroutine accelerating_pose(state, value, args)
    type(variational_state), intent(in) :: state
    real(real64), intent(out) :: value(:)
    class(*), intent(inout), optional :: args
    real(real64) :: rotation(3,3)

    rotation = state%orientation(1)%to_matrix()
    value(1) = state%position(1,1) - state%time**2
    value(2) = atan2(rotation(2,1), rotation(1,1)) - state%time**2
end subroutine

! ------------------------------------------------------------------------------
function test_variational_force_evaluation() result(rst)
    !! Verifies the three applied-load evaluation points with a scalar linear
    !! viscous damper.
    logical :: rst
    integer(int32), parameter :: mode(3) = [VI_FORCE_LEFT_ENDPOINT, &
        VI_FORCE_IMPLICIT_ENDPOINT, VI_FORCE_MIDPOINT]
    real(real64), parameter :: expected(3) = [0.8d0, 1.0d0 / 1.2d0, &
        0.9d0 / 1.1d0]
    type(rigid_body) :: bodies(1)
    type(variational_state) :: state
    type(variational_integrator) :: integrator
    procedure(variational_force), pointer :: force_ptr
    integer(int32) :: i

    rst = .true.
    bodies(1) = rigid_body(2.0d0)
    force_ptr => linear_damping_force
    do i = 1, size(mode)
        call initialize_variational_state(state, 1)
        state%velocity(1,1) = 1.0d0
        integrator%settings%force_evaluation = mode(i)
        call integrator%step(bodies, state, 0.1d0, &
            force_function = force_ptr)
        if (abs(state%velocity(1,1) - expected(i)) > 1.0d-9) then
            rst = .false.
            print "(A,I0)", &
                "TEST FAILED: test_variational_force_evaluation mode ", mode(i)
        end if
    end do
end function

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
    type(variational_integrator_info) :: info
        !! Convergence diagnostics.
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
    call integrator%step(bodies, state, dt, info = info)
    if (.not.info%converged .or. info%iterations /= 0 .or. &
        info%jacobian_singular) rst = .false.
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
function test_variational_recoverable_failure() result(rst)
    !! Verifies that requested diagnostics turn convergence failure into a
    !! recoverable return and preserve the last accepted state.
    logical :: rst
    type(rigid_body) :: bodies(1)
    type(variational_state) :: state
    type(variational_state), allocatable, dimension(:) :: solution
    type(variational_integrator) :: integrator
    type(variational_integrator_info) :: info
    real(real64), allocatable, dimension(:,:) :: multipliers
    procedure(variational_force), pointer :: force_ptr
    procedure(variational_constraint), pointer :: constraint_ptr
    procedure(variational_constraint_jacobian), pointer :: jacobian_ptr

    rst = .true.
    bodies(1) = rigid_body(2.0d0)
    call initialize_variational_state(state, 1)
    force_ptr => constant_force
    integrator%settings%maximum_iterations = 1
    integrator%settings%tolerance = 1.0d-30
    call integrator%step(bodies, state, 0.01d0, &
        force_function = force_ptr, info = info)
    if (info%converged .or. info%iterations /= 1 .or. &
        info%jacobian_singular .or. state%time /= 0.0d0 .or. &
        any(state%position /= 0.0d0) .or. any(state%velocity /= 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_variational_recoverable_failure - step"
    end if

    solution = integrator%solve(bodies, state, 0.01d0, 3, &
        force_function = force_ptr, info = info)
    if (info%converged .or. info%iterations /= 1 .or. &
        info%jacobian_singular .or. size(solution) /= 1 .or. &
        state%time /= 0.0d0) then
        rst = .false.
        print "(A)", "TEST FAILED: test_variational_recoverable_failure - solve"
    end if

    call initialize_variational_state(state, 1)
    constraint_ptr => late_inconsistent_constraint
    jacobian_ptr => zero_constraint_jacobian
    solution = integrator%solve(bodies, state, 0.01d0, 3, &
        constraint_count = 1, constraint = constraint_ptr, &
        constraint_jacobian = jacobian_ptr, multipliers = multipliers, &
        info = info)
    if (info%converged .or. info%iterations /= 1 .or. &
        .not.info%jacobian_singular .or. size(solution) /= 2 .or. &
        abs(state%time - 0.01d0) > epsilon(1.0d0) .or. &
        size(multipliers,1) /= 1 .or. size(multipliers,2) /= 1) then
        rst = .false.
        print "(A)", "TEST FAILED: test_variational_recoverable_failure - partial"
    end if
end function

! ------------------------------------------------------------------------------
subroutine late_inconsistent_constraint(state, value, args)
    type(variational_state), intent(in) :: state
    real(real64), intent(out), dimension(:) :: value
    class(*), intent(inout), optional :: args

    value = merge(1.0d0, 0.0d0, state%time > 0.015d0)
end subroutine

! ------------------------------------------------------------------------------
subroutine zero_constraint_jacobian(state, jacobian, args)
    type(variational_state), intent(in) :: state
    real(real64), intent(out), dimension(:,:) :: jacobian
    class(*), intent(inout), optional :: args

    jacobian = 0.0d0
end subroutine

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
    state%velocity(:,1) = 0.0d0
    state%angular_velocity(:,1) = 0.0d0
    integrator%settings%constraint_translation_scale = 1.0d3
    integrator%settings%constraint_rotation_scale = 0.25d0
    constraint_ptr => scaled_fixed_pose
    call integrator%step(bodies, state, 1.0d-3, 6, &
        constraint = constraint_ptr, multipliers = multipliers)
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
subroutine linear_damping_force(t, state, force, torque, args)
    real(real64), intent(in) :: t
    type(variational_state), intent(in) :: state
    real(real64), intent(out), dimension(:,:) :: force, torque
    class(*), intent(inout), optional :: args

    force = 0.0d0
    torque = 0.0d0
    force(1,1) = -4.0d0 * state%velocity(1,1)
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
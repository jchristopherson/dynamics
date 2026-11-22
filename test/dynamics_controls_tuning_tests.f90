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
module dynamics_controls_tuning_tests
    use iso_fortran_env
    use dynamics
    implicit none
contains
! ------------------------------------------------------------------------------
subroutine cl_output(K, Ki, C, D, r, x, xi, u, y)
    !! Evaluates the control input and plant output for a single-input,
    !! single-output plant.
    real(real64), intent(in), dimension(:,:) :: K, Ki, C, D
    real(real64), intent(in) :: r
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(in) :: xi
    real(real64), intent(out) :: u, y

    u = -dot_product(K(1,:), x) - Ki(1,1) * xi
    y = dot_product(C(1,:), x) + D(1,1) * u
end subroutine

! ------------------------------------------------------------------------------
subroutine cl_rhs(A, B, C, D, K, Ki, r, z, dz, e)
    !! The continuous-time closed-loop right-hand side for z = [x; xi].
    real(real64), intent(in), dimension(:,:) :: A, B, C, D, K, Ki
    real(real64), intent(in) :: r
    real(real64), intent(in), dimension(:) :: z
    real(real64), intent(out), dimension(:) :: dz
    real(real64), intent(out) :: e

    integer :: n
    real(real64) :: u, y

    n = size(A, 1)
    call cl_output(K, Ki, C, D, r, z(1:n), z(n + 1), u, y)
    dz(1:n) = matmul(A, z(1:n)) + B(:,1) * u
    dz(n + 1) = y - r
    e = y - r
end subroutine

! ------------------------------------------------------------------------------
function simulate_continuous_ise(A, B, C, D, K, Ki, r, tend, dt) result(ise)
    !! Integrates the nonlinear closed loop with RK4 and accumulates the
    !! integral of the squared error.
    real(real64), intent(in), dimension(:,:) :: A, B, C, D, K, Ki
    real(real64), intent(in) :: r, tend, dt
    real(real64) :: ise

    integer :: i, n, nsteps
    real(real64) :: e
    real(real64), allocatable, dimension(:) :: z, k1, k2, k3, k4, zt

    n = size(A, 1)
    allocate(z(n + 1), k1(n + 1), k2(n + 1), k3(n + 1), k4(n + 1), zt(n + 1))
    z = 0.0d0
    ise = 0.0d0
    nsteps = nint(tend / dt)
    do i = 1, nsteps
        call cl_rhs(A, B, C, D, K, Ki, r, z, k1, e)
        ise = ise + e * e * dt
        zt = z + 0.5d0 * dt * k1
        call cl_rhs(A, B, C, D, K, Ki, r, zt, k2, e)
        zt = z + 0.5d0 * dt * k2
        call cl_rhs(A, B, C, D, K, Ki, r, zt, k3, e)
        zt = z + dt * k3
        call cl_rhs(A, B, C, D, K, Ki, r, zt, k4, e)
        z = z + dt / 6.0d0 * (k1 + 2.0d0 * k2 + 2.0d0 * k3 + k4)
    end do
end function

! ------------------------------------------------------------------------------
function simulate_discrete_ise(A, B, C, D, K, Ki, Ts, r, nsteps) result(ise)
    !! Simulates the discrete closed loop and accumulates Ts * sum(e^2).
    real(real64), intent(in), dimension(:,:) :: A, B, C, D, K, Ki
    real(real64), intent(in) :: Ts, r
    integer, intent(in) :: nsteps
    real(real64) :: ise

    integer :: step
    real(real64) :: xi, u, y
    real(real64), allocatable, dimension(:) :: x

    allocate(x(size(A, 1)))
    x = 0.0d0
    xi = 0.0d0
    ise = 0.0d0
    do step = 1, nsteps
        call cl_output(K, Ki, C, D, r, x, xi, u, y)
        ise = ise + Ts * (y - r)**2
        x = matmul(A, x) + B(:,1) * u
        xi = xi + Ts * (y - r)
    end do
end function

! ------------------------------------------------------------------------------
subroutine build_test_plant(sys_c, sys_d, K, Ki, ctrl)
    !! Builds a two-state plant with a feedthrough term (both continuous and
    !! discrete) and a stable starting controller.
    type(state_space), intent(out) :: sys_c
    type(discrete_state_space), intent(out) :: sys_d
    real(real64), intent(out), allocatable, dimension(:,:) :: K, Ki
    type(integral_state_feedback_controller), intent(out) :: ctrl

    real(real64), allocatable, dimension(:,:) :: A, B, C, D
    allocate(A(2,2), B(2,1), C(1,2), D(1,1))
    A = reshape([0.0d0, -2.0d0, 1.0d0, -3.0d0], [2, 2])
    B = reshape([0.0d0, 1.0d0], [2, 1])
    C = reshape([1.0d0, 0.0d0], [1, 2])
    D = reshape([0.2d0], [1, 1])
    sys_c = state_space(A, B, C, D)
    K = reshape([2.0d0, 1.0d0], [1, 2])
    Ki = reshape([1.0d0], [1, 1])
    ctrl = integral_state_feedback_controller(K, Ki)
end subroutine

! ------------------------------------------------------------------------------
function test_isf_continuous_ise() result(rst)
    !! The closed-form ISE must agree with a direct simulation of the closed
    !! loop, including a nonzero feedthrough term.
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-3

    ! Local Variables
    real(real64) :: ise_cf, ise_sim
    real(real64), allocatable, dimension(:,:) :: K, Ki, A, B, C, D
    type(state_space) :: sys_c
    type(discrete_state_space) :: sys_d
    type(integral_state_feedback_controller) :: ctrl

    call build_test_plant(sys_c, sys_d, K, Ki, ctrl)
    A = sys_c%A
    B = sys_c%B
    C = sys_c%C
    D = sys_c%D

    rst = ctrl%is_stable(sys_c)
    ise_cf = ctrl%integral_squared_error(sys_c)
    ise_sim = simulate_continuous_ise(A, B, C, D, K, Ki, 1.0d0, 30.0d0, 1.0d-4)
    rst = rst .and. abs(ise_cf - ise_sim) <= tol * ise_sim
end function

! ------------------------------------------------------------------------------
function test_isf_discrete_ise() result(rst)
    !! The closed-form discrete ISE must agree with the summed squared error
    !! from a direct simulation of the discrete closed loop.
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    real(real64) :: ise_cf, ise_sim, Ts
    real(real64), allocatable, dimension(:,:) :: K, Ki, A, B, C, D
    type(discrete_state_space) :: sys_d
    type(integral_state_feedback_controller) :: ctrl

    A = reshape([0.9d0, 0.0d0, 0.1d0, 0.8d0], [2, 2])
    B = reshape([0.0d0, 0.2d0], [2, 1])
    C = reshape([1.0d0, 0.0d0], [1, 2])
    D = reshape([0.1d0], [1, 1])
    Ts = 0.1d0
    sys_d = discrete_state_space(A, B, C, D, Ts)
    K = reshape([1.0d0, 0.5d0], [1, 2])
    Ki = reshape([0.3d0], [1, 1])
    ctrl = integral_state_feedback_controller(K, Ki)

    rst = ctrl%is_stable(sys_d)
    ise_cf = ctrl%integral_squared_error(sys_d)
    ise_sim = simulate_discrete_ise(A, B, C, D, K, Ki, Ts, 1.0d0, 5000)
    rst = rst .and. abs(ise_cf - ise_sim) <= tol * ise_sim
end function

! ------------------------------------------------------------------------------
function test_isf_unstable_gains() result(rst)
    !! Unstable closed loops must be reported as unstable with an infinite ISE.
    logical :: rst

    ! Local Variables
    real(real64) :: ise
    real(real64), allocatable, dimension(:,:) :: K, Ki
    type(state_space) :: sys_c
    type(discrete_state_space) :: sys_d
    type(integral_state_feedback_controller) :: ctrl

    call build_test_plant(sys_c, sys_d, K, Ki, ctrl)
    K = reshape([-50.0d0, 0.0d0], [1, 2])
    ctrl = integral_state_feedback_controller(K, Ki)

    ise = ctrl%integral_squared_error(sys_c)
    rst = .not.ctrl%is_stable(sys_c) .and. ise > huge(ise)
end function

! ------------------------------------------------------------------------------
function test_isf_optimizer_continuous() result(rst)
    !! Optimizing the continuous-time gains must keep the loop stable, reduce
    !! the ISE, and report the ISE at the optimized gains.
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    real(real64) :: ise0, ise1, fout
    real(real64), allocatable, dimension(:,:) :: K, Ki
    type(state_space) :: sys_c
    type(discrete_state_space) :: sys_d
    type(integral_state_feedback_controller) :: ctrl
    type(controller_optimizer) :: opt

    call build_test_plant(sys_c, sys_d, K, Ki, ctrl)
    ise0 = ctrl%integral_squared_error(sys_c)

    call opt%optimize(sys_c, ctrl, fout)
    ise1 = ctrl%integral_squared_error(sys_c)
    rst = ctrl%is_stable(sys_c) .and. ise1 < ise0
    rst = rst .and. abs(ise1 - fout) <= tol * max(1.0d0, fout)
end function

! ------------------------------------------------------------------------------
function test_isf_optimizer_discrete() result(rst)
    !! Optimizing the discrete-time gains must keep the loop stable and reduce
    !! the ISE.
    logical :: rst

    ! Local Variables
    real(real64) :: ise0, ise1, fout
    real(real64), allocatable, dimension(:,:) :: K, Ki, A, B, C, D
    type(discrete_state_space) :: sys_d
    type(integral_state_feedback_controller) :: ctrl
    type(controller_optimizer) :: opt

    A = reshape([0.9d0, 0.0d0, 0.1d0, 0.8d0], [2, 2])
    B = reshape([0.0d0, 0.2d0], [2, 1])
    C = reshape([1.0d0, 0.0d0], [1, 2])
    D = reshape([0.1d0], [1, 1])
    sys_d = discrete_state_space(A, B, C, D, 0.1d0)
    K = reshape([1.0d0, 0.5d0], [1, 2])
    Ki = reshape([0.3d0], [1, 1])
    ctrl = integral_state_feedback_controller(K, Ki)
    ise0 = ctrl%integral_squared_error(sys_d)

    call opt%optimize(sys_d, ctrl, fout)
    ise1 = ctrl%integral_squared_error(sys_d)
    rst = ctrl%is_stable(sys_d) .and. ise1 < ise0
    rst = rst .and. abs(ise1 - fout) <= 1.0d-6 * max(1.0d0, fout)
end function

! ------------------------------------------------------------------------------
function test_isf_optimizer_control_weight() result(rst)
    !! A control rate weight must reduce the magnitude of the optimized gains
    !! relative to the ISE-only optimum, and keep the loop stable.
    logical :: rst

    ! Local Variables
    real(real64) :: fout
    real(real64), allocatable, dimension(:,:) :: K, Ki
    type(state_space) :: sys_c
    type(discrete_state_space) :: sys_d
    type(integral_state_feedback_controller) :: ctrl_ise, ctrl_pen
    type(controller_optimizer) :: opt

    call build_test_plant(sys_c, sys_d, K, Ki, ctrl_ise)
    ctrl_pen = ctrl_ise

    call opt%optimize(sys_c, ctrl_ise, fout)
    opt%control_weight = 10.0d0
    call opt%optimize(sys_c, ctrl_pen, fout)

    rst = ctrl_pen%is_stable(sys_c)
    rst = rst .and. sum(ctrl_pen%K**2) + sum(ctrl_pen%Ki**2) < &
        sum(ctrl_ise%K**2) + sum(ctrl_ise%Ki**2)
end function

! ------------------------------------------------------------------------------
function test_isf_apply_to_ss_plant() result(rst)
    !! The closed loop realized by apply_to_state_space_plant must reproduce
    !! the plant and controller equations for a two-output plant with
    !! feedthrough.  The inputs are the references and the outputs are
    !! [u; y].
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-12

    ! Local Variables
    real(real64) :: u
    real(real64), allocatable, dimension(:,:) :: A, B, C, D, K, Ki
    real(real64), allocatable, dimension(:) :: x, xi, z, r, y, dx, dxi, dz, out, &
        expected_dz, expected_out
    type(state_space) :: sys_c, sys_r
    type(integral_state_feedback_controller) :: ctrl

    A = reshape([0.0d0, 1.0d0, -2.0d0, -3.0d0], [2, 2])
    B = reshape([0.0d0, 1.0d0], [2, 1])
    C = reshape([1.0d0, 0.0d0, 0.0d0, 1.0d0], [2, 2])
    D = reshape([0.2d0, 0.1d0], [2, 1])
    sys_c = state_space(A, B, C, D)
    K = reshape([2.0d0, 1.0d0], [1, 2])
    Ki = reshape([1.0d0, 0.5d0], [1, 2])
    ctrl = integral_state_feedback_controller(K, Ki)
    sys_r = ctrl%apply_to_state_space_plant(sys_c)

    rst = .false.
    if (size(sys_r%A, 1) /= 4 .or. size(sys_r%B, 2) /= 2 .or. &
        size(sys_r%C, 1) /= 3 .or. size(sys_r%D, 2) /= 2) return

    z = [0.3d0, -0.7d0, 0.4d0, -0.2d0]
    r = [0.5d0, -0.1d0]
    x = z(1:2)
    xi = z(3:4)
    u = -dot_product(K(1,:), x) - dot_product(Ki(1,:), xi)
    y = matmul(C, x) + D(:,1) * u
    dx = matmul(A, x) + B(:,1) * u
    dxi = y - r
    expected_dz = [dx, dxi]
    expected_out = [u, y]

    dz = matmul(sys_r%A, z) + matmul(sys_r%B, r)
    out = matmul(sys_r%C, z) + matmul(sys_r%D, r)
    rst = all(abs(dz - expected_dz) <= tol) .and. all(abs(out - expected_out) <= tol)
end function

! ------------------------------------------------------------------------------
function test_isf_apply_to_dss_plant() result(rst)
    !! The discrete closed loop realized by apply_to_discrete_state_space_plant
    !! must reproduce the discrete plant and integrator updates for a
    !! two-output plant with feedthrough, and keep the sample time.
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-12
    real(real64), parameter :: Ts = 0.1d0

    ! Local Variables
    real(real64) :: u
    real(real64), allocatable, dimension(:,:) :: A, B, C, D, K, Ki
    real(real64), allocatable, dimension(:) :: x, xi, z, r, y, dx, dxi, dz, out, &
        expected_dz, expected_out
    type(discrete_state_space) :: sys_d, sys_r
    type(integral_state_feedback_controller) :: ctrl

    A = reshape([0.9d0, 0.0d0, 0.1d0, 0.8d0], [2, 2])
    B = reshape([0.0d0, 0.2d0], [2, 1])
    C = reshape([1.0d0, 0.0d0, 0.0d0, 1.0d0], [2, 2])
    D = reshape([0.1d0, 0.05d0], [2, 1])
    sys_d = discrete_state_space(A, B, C, D, Ts)
    K = reshape([1.0d0, 0.5d0], [1, 2])
    Ki = reshape([0.3d0, 0.2d0], [1, 2])
    ctrl = integral_state_feedback_controller(K, Ki)
    sys_r = ctrl%apply_to_discrete_state_space_plant(sys_d)

    rst = .false.
    if (size(sys_r%A, 1) /= 4 .or. size(sys_r%B, 2) /= 2 .or. &
        size(sys_r%C, 1) /= 3 .or. size(sys_r%D, 2) /= 2) return
    if (sys_r%sample_time /= Ts) return

    z = [0.3d0, -0.7d0, 0.4d0, -0.2d0]
    r = [0.5d0, -0.1d0]
    x = z(1:2)
    xi = z(3:4)
    u = -dot_product(K(1,:), x) - dot_product(Ki(1,:), xi)
    y = matmul(C, x) + D(:,1) * u
    dx = matmul(A, x) + B(:,1) * u
    dxi = xi + Ts * (y - r)
    expected_dz = [dx, dxi]
    expected_out = [u, y]

    dz = matmul(sys_r%A, z) + matmul(sys_r%B, r)
    out = matmul(sys_r%C, z) + matmul(sys_r%D, r)
    rst = all(abs(dz - expected_dz) <= tol) .and. all(abs(out - expected_out) <= tol)
end function

end module

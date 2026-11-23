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
module dynamics_state_space_tests
    use iso_fortran_env
    use fortran_test_helper
    use diffeq, only : runge_kutta_45, ode_container
    use dynamics
    implicit none

    real(real64), private, parameter :: m = 6.0d-4
    real(real64), private, parameter :: k = 2.0d6
    real(real64), private, parameter :: zeta = 1.0d-1
    real(real64), private, parameter :: b = 2.0d0 * zeta * sqrt(k * m)
    real(real64), private, parameter :: U = 1.0d4
    real(real64), private, parameter :: w = 5.0d2

contains
! ------------------------------------------------------------------------------
subroutine ode_fcn(t, x, dxdt, args)
    real(real64), intent(in) :: t, x(:)
    real(real64), intent(out) :: dxdt(:)
    class(*), intent(inout), optional :: args

    ! m x" + b x' + k x = U sin(w * t)
    dxdt(1) = x(2)
    dxdt(2) = (1.0d0 / m) * (U * sin(w * t) - b * x(2) - k * x(1))
end subroutine

subroutine ode_excite_fcn(t, f, args)
    real(real64), intent(in) :: t
    real(real64), intent(out) :: f(:)
    class(*), intent(in), optional :: args
    f = U * sin(w * t)
end subroutine

function test_lti_solve() result(rst)
    logical :: rst

    ! Parameters
    integer, parameter :: ntime = 1000
    real(real64), parameter :: tmax = 1.0d-1
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    integer :: i
    real(real64) :: dt, t(ntime), ic(2)
    real(real64), allocatable, dimension(:,:) :: sol, ltiSol
    procedure(ss_excitation), pointer :: excite
    type(ode_container) :: odeMdl
    type(state_space) :: mdl
    type(runge_kutta_45) :: integrator

    ! Initialization
    rst = .true.
    excite => ode_excite_fcn
    odeMdl%fcn => ode_fcn

    ! Build the state space model
    mdl = state_space(m, b, k)
    mdl%C(1, 2) = 0.0d0

    ! Build a time vector and an initial conditions vector
    dt = tmax / (ntime - 1.0d0)
    t = (/ (i * dt, i = 0, ntime - 1) /)
    ic = 0.0d0

    ! Solve using LTI
    ltiSol = lti_solve(mdl, excite, t, ic)

    ! Solve using a traditional ODE solver
    call integrator%solve(odeMdl, t, ic)
    sol = integrator%get_solution()

    ! Test
    if (.not.assert(ltiSol, sol(:,1:2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_lti_solve -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_state_space_initialize() result(rst)
    logical :: rst

    ! Local Variables
    real(real64) :: m1, m2, b1, b2, k1, k2, M(2,2), B(2,2), K(2,2), Aans(4,4), &
        Bans(4,2), Cans(1,4), Dans(1,2)
    type(state_space) :: mdl
    
    ! Initialization
    rst = .true.
    call random_number(m1)
    call random_number(m2)
    call random_number(b1)
    call random_number(b2)
    call random_number(k1)
    call random_number(k2)

    M = reshape([ &
        m1, 0.0d0, &
        0.0d0, m2 &
        ], &
    [2, 2])
    B = reshape([ &
        b1+b2, -b2, &
        -b2, b2 &
        ], &
    [2, 2])
    K = reshape([ &
        k1+k2, -k2, &
        -k2, k2 &
        ], &
    [2, 2])

    Aans = reshape([ &
        0.0d0, 0.0d0, -(k1+k2)/m1, k2/m2, &
        0.0d0, 0.0d0, k2/m1, -k2/m2, &
        1.0d0, 0.0d0, -(b1+b2)/m1, b2/m2, &
        0.0d0, 1.0d0, b2/m1, -b2/m2 &
        ], &
    [4, 4])
    Bans = reshape([ &
        0.0d0, 0.0d0, 1.0d0/m1, 0.0d0, &
        0.0d0, 0.0d0, 0.0d0, 1.0d0/m2 &
        ], &
    [4, 2])
    Cans = 1.0d0
    Dans = 0.0d0

    ! Create a state space model
    mdl = state_space(M, B, K)

    ! Tests
    if (.not.assert(mdl%A, Aans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_state_space_initialize -1"
    end if
    if (.not.assert(mdl%B, Bans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_state_space_initialize -2"
    end if
    if (.not.assert(mdl%C, Cans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_state_space_initialize -3"
    end if
    if (.not.assert(mdl%D, Dans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_state_space_initialize -4"
    end if
end function

! ------------------------------------------------------------------------------
function test_state_space_poles_zeros() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    complex(real64), allocatable, dimension(:) :: pAns, zAns, poles, zeros
    real(real64) :: m, b, k
    type(transfer_function) :: tf
    type(state_space) :: mdl

    ! Initialization
    rst = .true.
    call random_number(m)
    call random_number(b)
    call random_number(k)
    tf = transfer_function([k, b], [k, b, m])
    mdl = tf%to_ccf_state_space()

    ! Compute the solution
    pAns = tf%poles()
    zAns = tf%zeros()

    ! Test
    poles = mdl%poles()
    if (.not.assert(poles, pAns, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_state_space_poles_zeros -1"
        print *, "ANSWER: ", pAns
        print *, "COMPUTED: ", poles
    end if

    zeros = mdl%zeros()
    if (.not.assert(zeros, zAns, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_state_space_poles_zeros -2"
        print *, "ANSWER: ", zAns
        print *, "COMPUTED: ", zeros
    end if
end function

! ------------------------------------------------------------------------------
function test_state_space_pid_order() result(rst)
    logical :: rst

    integer, parameter :: n = 3
    real(real64), parameter :: kp = 2.0d0, ki = 3.0d0
    real(real64), parameter :: kd = 4.0d0, tau = 2.0d0
    real(real64) :: a(n,n), b(n,1), c(1,n), d(1,1)
    real(real64) :: expected_a(n+2,n+2), expected_b(n+2,1)
    real(real64) :: expected_c(1,n+2), expected_d(1,1)
    type(state_space) :: mdl

    rst = .true.
    a = 0.0d0
    b(:,1) = [1.0d0, 2.0d0, 3.0d0]
    c(1,:) = [4.0d0, 5.0d0, 6.0d0]
    d = 0.0d0
    expected_a = 0.0d0
    expected_b = 0.0d0
    expected_c = 0.0d0
    expected_d = 0.0d0
    expected_a(1:n,1:n) = a - (kp + kd / tau) * matmul(b,c)
    expected_a(n+1,1:n) = -c(1,:)
    expected_a(n+2,1:n) = expected_a(n+1,1:n) / tau
    expected_a(1:n,n+1) = ki * b(:,1)
    expected_a(1:n,n+2) = -(kd / tau) * b(:,1)
    expected_a(n+2,n+2) = -1.0d0 / tau
    expected_b(1:n,1) = (kp + kd / tau) * b(:,1)
    expected_b(n+1,1) = 1.0d0
    expected_b(n+2,1) = 1.0d0 / tau
    expected_c(1,1:n) = c(1,:)

    mdl = state_space(kp, ki, kd, tau, a, b, c, d)
    if (maxval(abs(mdl%A - expected_a)) > 1.0d-12) rst = .false.
    if (maxval(abs(mdl%B - expected_b)) > 1.0d-12) rst = .false.
    if (maxval(abs(mdl%C - expected_c)) > 1.0d-12) rst = .false.
    if (maxval(abs(mdl%D - expected_d)) > 1.0d-12) rst = .false.
    if (.not.rst) print "(A)", &
        "TEST FAILED: test_state_space_pid_order"
end function

! ------------------------------------------------------------------------------
function test_state_space_to_transfer_function() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-8
    complex(real64), parameter :: j = (0.0d0, 1.0d0)

    ! Local Variables
    real(real64) :: m, b, k, omega
    complex(real64) :: s, ans
    complex(real64), allocatable, dimension(:,:) :: z
    type(transfer_function) :: tf
    type(state_space) :: mdl

    ! Initialization
    rst = .true.
    call random_number(m)
    call random_number(b)
    call random_number(k)
    call random_number(omega)
    s = j * omega
    tf = transfer_function([k, b], [k, b, m])
    mdl = tf%to_ccf_state_space()

    ! Test
    ans = tf%evaluate(s)
    z = mdl%transfer_function(s)
    if (.not.assert(ans, z(1,1))) then
        rst = .false.
    end if
end function

! ------------------------------------------------------------------------------
function test_state_space_to_transfer_function_round_trip() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-8
    complex(real64), parameter :: s = (0.3d0, 1.7d0)

    ! Local Variables
    complex(real64) :: ans, z
    type(transfer_function) :: tf, tf2
    type(state_space) :: mdl

    ! Initialization
    rst = .true.
    tf = transfer_function([4.0d0, 3.0d0], [4.0d0, 3.0d0, 2.0d0])
    mdl = tf%to_ccf_state_space()

    ! Process
    tf2 = mdl%to_transfer_function()

    ! Test
    ans = tf%evaluate(s)
    z = tf2%evaluate(s)
    if (.not.assert(ans, z, tol)) rst = .false.
end function

! ------------------------------------------------------------------------------
function test_matrix_exponential() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-12

    ! Local Variables
    real(real64), allocatable, dimension(:,:) :: a, ans, z

    ! Initialization
    rst = .true.

    ! Test a nilpotent matrix: exp([0 2; 0 0]) = [1 2; 0 1]
    a = reshape([0.0d0, 0.0d0, 2.0d0, 0.0d0], [2, 2])
    ans = reshape([1.0d0, 0.0d0, 2.0d0, 1.0d0], [2, 2])
    z = matrix_exponential(a)
    if (.not.assert(z, ans, tol)) rst = .false.

    ! Test a diagonal matrix
    a = reshape([-1.0d0, 0.0d0, 0.0d0, -2.0d0], [2, 2])
    ans = reshape([exp(-1.0d0), 0.0d0, 0.0d0, exp(-2.0d0)], [2, 2])
    z = matrix_exponential(a)
    if (.not.assert(z, ans, tol)) rst = .false.
end function

! ------------------------------------------------------------------------------
function test_discrete_state_space_to_transfer_function() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-12

    ! Local Variables
    complex(real64), allocatable, dimension(:) :: p, z
    type(discrete_state_space) :: mdl
    type(discrete_transfer_function) :: tf

    ! Initialization
    rst = .true.
    mdl = discrete_state_space(reshape([0.5d0], [1, 1]), &
        reshape([1.0d0], [1, 1]), reshape([1.0d0], [1, 1]), &
        reshape([2.0d0], [1, 1]), 0.1d0)

    ! Process
    tf = mdl%to_transfer_function()
    p = tf%poles()
    z = tf%zeros()

    ! Test: H(z) = 2z / (z - 0.5)
    if (size(p) /= 1 .or. size(z) /= 1) then
        rst = .false.
        return
    end if
    if (.not.assert(real(p(1)), 0.5d0, tol)) rst = .false.
    if (.not.assert(abs(z(1)), 0.0d0, tol)) rst = .false.
end function

! ------------------------------------------------------------------------------
function test_discrete_state_space_zeros_evaluate() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: ts = 0.1d0
    real(real64), parameter :: tol = 1.0d-10
    real(real64), parameter :: tol_eval = 1.0d-12

    ! Local Variables
    complex(real64) :: z
    complex(real64), allocatable, dimension(:) :: zeros
    real(real64), allocatable, dimension(:) :: x_next, y
    integer :: i
    type(discrete_state_space) :: mdl

    ! Initialization
    rst = .true.

    ! Test: H(z) = 2z / (z - 0.5) has a single zero at the origin
    mdl = discrete_state_space(reshape([0.5d0], [1, 1]), &
        reshape([1.0d0], [1, 1]), reshape([1.0d0], [1, 1]), &
        reshape([2.0d0], [1, 1]), ts)
    zeros = mdl%zeros()
    if (size(zeros) /= 1) then
        rst = .false.
    else if (abs(zeros(1)) > tol) then
        rst = .false.
    end if
    ! x_{k+1} = 0.5(0.3) + 1.5 = 1.65 and y_k = 0.3 + 2(1.5) = 3.3
    x_next = mdl%evaluate_next_state([1.5d0], [0.3d0])
    y = mdl%evaluate_output([1.5d0], [0.3d0])
    if (abs(x_next(1) - 1.65d0) > tol_eval) rst = .false.
    if (abs(y(1) - 3.3d0) > tol_eval) rst = .false.

    ! Test: H(z) = (0.5z^2 + 2.85z - 0.85) / ((z - 0.5)(z + 0.2)) + 0.5,
    ! which has two zeros
    mdl = discrete_state_space(reshape([0.5d0, 0.0d0, 0.0d0, -0.2d0], [2, 2]), &
        reshape([1.0d0, 1.0d0], [2, 1]), reshape([1.0d0, 2.0d0], [1, 2]), &
        reshape([0.5d0], [1, 1]), ts)
    zeros = mdl%zeros()
    if (size(zeros) /= 2) then
        rst = .false.
    else
        do i = 1, size(zeros)
            z = zeros(i)
            if (abs(0.5d0 * z * z + 2.85d0 * z - 0.85d0) > tol) rst = .false.
        end do
    end if
    ! x_{k+1} = [0.5(1) + 3; -0.2(2) + 3] = [3.5; 2.6]
    ! y_k = (1 + 2(2)) + 0.5(3) = 6.5
    x_next = mdl%evaluate_next_state([3.0d0], [1.0d0, 2.0d0])
    y = mdl%evaluate_output([3.0d0], [1.0d0, 2.0d0])
    if (abs(x_next(1) - 3.5d0) > tol_eval) rst = .false.
    if (abs(x_next(2) - 2.6d0) > tol_eval) rst = .false.
    if (abs(y(1) - 6.5d0) > tol_eval) rst = .false.

    ! Test: H(z) = 1 / (z - 0.5) has no finite zeros
    mdl = discrete_state_space(reshape([0.5d0], [1, 1]), &
        reshape([1.0d0], [1, 1]), reshape([1.0d0], [1, 1]), &
        reshape([0.0d0], [1, 1]), ts)
    zeros = mdl%zeros()
    if (size(zeros) /= 0) rst = .false.

    if (.not.rst) print "(A)", "TEST FAILED: test_discrete_state_space_zeros_evaluate"
end function

! ------------------------------------------------------------------------------
function test_discrete_canonical_forms() result(rst)
    logical :: rst
    real(real64), parameter :: tol = 1.0d-12
    real(real64), dimension(2,2) :: expected_a
    real(real64), dimension(2,1) :: expected_b
    real(real64), dimension(1,2) :: expected_c
    type(discrete_transfer_function) :: tf
    type(discrete_state_space) :: ccf
    type(discrete_state_space) :: ocf

    rst = .true.
    tf = discrete_transfer_function([1.0d0, 0.6d0], &
        [2.0d0, -1.2d0, 0.2d0], 0.125d0)
    ccf = tf%to_ccf_state_space()
    ocf = tf%to_ocf_state_space()
    expected_a = reshape([0.0d0, -0.1d0, 1.0d0, 0.6d0], [2,2])
    expected_b = reshape([0.0d0, 1.0d0], [2,1])
    expected_c = reshape([-0.05d0, 0.6d0], [1,2])
    if (.not.assert(ccf%A, expected_a, tol)) rst = .false.
    if (.not.assert(ccf%B, expected_b, tol)) rst = .false.
    if (.not.assert(ccf%C, expected_c, tol)) rst = .false.
    if (abs(ccf%D(1,1) - 0.5d0) > tol) rst = .false.
    if (.not.assert(ocf%A, transpose(expected_a), tol)) rst = .false.
    if (.not.assert(ocf%B, transpose(expected_c), tol)) rst = .false.
    if (.not.assert(ocf%C, transpose(expected_b), tol)) rst = .false.
    if (abs(ocf%D(1,1) - 0.5d0) > tol) rst = .false.

    if (.not.check_discrete_canonical([1.0d0, 0.6d0], &
        [2.0d0, -1.2d0, 0.2d0], 2)) rst = .false.
    if (.not.check_discrete_canonical([0.0d0, 1.0d0, -0.5d0], &
        [2.0d0], 2)) rst = .false.
    if (.not.check_discrete_canonical([1.0d0, 0.2d0, -0.3d0, 0.4d0], &
        [2.0d0, -0.4d0], 3)) rst = .false.
    if (.not.check_discrete_canonical([0.0d0, 0.3d0], &
        [1.0d0, -1.2d0, 0.36d0], 2)) rst = .false.
    if (.not.check_discrete_canonical([3.0d0], &
        [2.0d0, -1.0d0, 0.0d0, 0.0d0], 1)) rst = .false.
    if (.not.check_discrete_canonical([0.0d0, 0.0d0], [2.0d0], 0)) rst = .false.
    if (.not.check_discrete_canonical([3.0d0], [2.0d0], 0)) rst = .false.
    if (.not.rst) print "(A)", "TEST FAILED: test_discrete_canonical_forms"
end function

! ------------------------------------------------------------------------------
function check_discrete_canonical(b, a, n) result(rst)
    real(real64), intent(in), dimension(:) :: b
    real(real64), intent(in), dimension(:) :: a
    integer, intent(in) :: n
    logical :: rst
    integer, parameter :: ns = 20
    real(real64), parameter :: ts = 0.125d0
    real(real64), parameter :: tol = 1.0d-11
    integer :: i
    integer :: j
    integer :: form
    real(real64) :: value
    real(real64), dimension(ns) :: u
    real(real64), dimension(ns) :: y
    real(real64), allocatable, dimension(:) :: x
    type(discrete_transfer_function) :: tf
    type(discrete_transfer_function) :: round_trip
    type(discrete_state_space) :: ss

    rst = .true.
    tf = discrete_transfer_function(b, a, ts)
    u = [(sin(real(i, real64)) + 0.2d0 * cos(0.7d0 * real(i, real64)), i = 1, ns)]
    y = 0.0d0
    do i = 1, ns
        do j = 1, min(i, size(b))
            y(i) = y(i) + b(j) * u(i - j + 1)
        end do
        do j = 2, min(i, size(a))
            y(i) = y(i) - a(j) * y(i - j + 1)
        end do
        y(i) = y(i) / a(1)
    end do
    do form = 1, 2
        if (form == 1) then
            ss = tf%to_ccf_state_space()
        else
            ss = tf%to_ocf_state_space()
        end if
        if (size(ss%A, 1) /= n .or. size(ss%A, 2) /= n .or. &
            size(ss%B, 1) /= n .or. size(ss%B, 2) /= 1 .or. &
            size(ss%C, 1) /= 1 .or. size(ss%C, 2) /= n .or. &
            size(ss%D, 1) /= 1 .or. size(ss%D, 2) /= 1) then
            rst = .false.
            return
        end if
        if (ss%sample_time /= ts) rst = .false.
        allocate(x(n), source = 0.0d0)
        do i = 1, ns
            value = dot_product(ss%C(1,:), x) + ss%D(1,1) * u(i)
            if (abs(value - y(i)) > tol) rst = .false.
            x = matmul(ss%A, x) + ss%B(:,1) * u(i)
        end do
        deallocate(x)
        round_trip = ss%to_transfer_function()
        do i = 1, n + 1
            value = 0.0d0
            if (i <= size(a)) value = a(i) / a(1)
            if (abs(round_trip%X%get(i) - value) > tol) rst = .false.
            value = 0.0d0
            if (i <= size(b)) value = b(i) / a(1)
            if (abs(round_trip%Y%get(i) - value) > tol) rst = .false.
        end do
    end do
end function

! ------------------------------------------------------------------------------
function test_discrete_tf_evaluate() result(rst)
    logical :: rst
    real(real64), parameter :: ts = 0.125d0
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: tol = 1.0d-12
    integer :: i
    real(real64), dimension(4) :: omega
    complex(real64) :: q
    complex(real64) :: expected
    complex(real64), dimension(4) :: actual
    type(discrete_transfer_function) :: tf

    rst = .true.
    omega = [0.0d0, 2.4d0, -3.1d0, pi / ts]
    tf = discrete_transfer_function([0.4d0, -0.2d0, 0.1d0], &
        [2.0d0, -0.8d0], ts)
    actual = tf%evaluate(omega)
    do i = 1, size(omega)
        q = cmplx(cos(omega(i) * ts), -sin(omega(i) * ts), kind = real64)
        expected = (0.4d0 - 0.2d0 * q + 0.1d0 * q**2) / (2.0d0 - 0.8d0 * q)
        if (abs(actual(i) - expected) > tol) rst = .false.
        if (abs(tf%evaluate(omega(i)) - expected) > tol) rst = .false.
    end do
    tf = discrete_transfer_function([0.0d0, 1.0d0], [1.0d0], ts)
    q = cmplx(cos(omega(2) * ts), -sin(omega(2) * ts), kind = real64)
    if (abs(tf%evaluate(omega(2)) - q) > tol) rst = .false.
    if (abs(tf%evaluate(omega(2) + 2.0d0 * pi / ts) - q) > tol) rst = .false.
    tf = discrete_transfer_function([3.0d0], [2.0d0], ts)
    if (abs(tf%evaluate(omega(3)) - 1.5d0) > tol) rst = .false.
    tf = discrete_transfer_function([0.3d0], [1.0d0, -0.4d0, 0.15d0], ts)
    q = cmplx(cos(omega(3) * ts), -sin(omega(3) * ts), kind = real64)
    expected = 0.3d0 / (1.0d0 - 0.4d0 * q + 0.15d0 * q**2)
    if (abs(tf%evaluate(omega(3)) - expected) > tol) rst = .false.
    if (.not.rst) print "(A)", "TEST FAILED: test_discrete_tf_evaluate"
end function

! ------------------------------------------------------------------------------
subroutine invalid_discrete_canonical(selection)
    character(len = *), intent(in) :: selection
    type(discrete_transfer_function) :: tf
    type(discrete_state_space) :: ss
    type(discrete_state_space) :: mdl
    real(real64), allocatable, dimension(:) :: y

    select case (selection)
    case ("discrete-invalid-ccf")
        tf = discrete_transfer_function([1.0d0], [0.0d0, 1.0d0], 0.1d0)
        ss = tf%to_ccf_state_space()
    case ("discrete-invalid-ocf")
        tf = discrete_transfer_function([1.0d0], [0.0d0, 1.0d0], 0.1d0)
        ss = tf%to_ocf_state_space()
    case ("discrete-invalid-sample-time")
        tf = discrete_transfer_function([1.0d0], [1.0d0], 0.1d0)
        tf%sample_time = 0.0d0
        ss = tf%to_ccf_state_space()
    case ("discrete-invalid-uninitialized")
        tf%sample_time = 0.1d0
        ss = tf%to_ocf_state_space()
    case ("discrete-invalid-evaluate-input")
        mdl = discrete_state_space(reshape([0.5d0, 0.0d0, 0.0d0, 0.5d0], [2, 2]), &
            reshape([1.0d0, 1.0d0], [2, 1]), reshape([1.0d0, 1.0d0], [1, 2]), &
            reshape([0.0d0], [1, 1]), 0.1d0)
        y = mdl%evaluate_next_state([1.0d0, 2.0d0], [1.0d0, 2.0d0])
    case ("discrete-invalid-evaluate-state")
        mdl = discrete_state_space(reshape([0.5d0, 0.0d0, 0.0d0, 0.5d0], [2, 2]), &
            reshape([1.0d0, 1.0d0], [2, 1]), reshape([1.0d0, 1.0d0], [1, 2]), &
            reshape([0.0d0], [1, 1]), 0.1d0)
        y = mdl%evaluate_output([1.0d0], [1.0d0])
    end select
end subroutine

! ------------------------------------------------------------------------------
function test_discretize_zoh_scalar() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-12

    ! Local Variables
    type(discrete_state_space) :: dsys
    type(state_space) :: sys

    ! Initialization
    rst = .true.
    sys = state_space(reshape([-1.0d0], [1, 1]), reshape([1.0d0], [1, 1]), &
        reshape([1.0d0], [1, 1]), reshape([0.0d0], [1, 1]))

    ! Process
    dsys = c2d(sys, 0.5d0)

    ! Test: Ad = exp(-T), Bd = 1 - exp(-T)
    if (.not.assert(dsys%A, reshape([exp(-0.5d0)], [1, 1]), tol)) rst = .false.
    if (.not.assert(dsys%B, reshape([1.0d0 - exp(-0.5d0)], [1, 1]), tol)) &
        rst = .false.
end function

! ------------------------------------------------------------------------------
function test_discretize_tustin_scalar() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-12

    ! Local Variables
    type(discrete_state_space) :: dsys
    type(state_space) :: sys

    ! Initialization
    rst = .true.
    sys = state_space(reshape([-1.0d0], [1, 1]), reshape([1.0d0], [1, 1]), &
        reshape([1.0d0], [1, 1]), reshape([0.0d0], [1, 1]))

    ! Process
    dsys = c2d(sys, 0.5d0, DISCRETIZE_TUSTIN)

    ! Test: Ad = 0.6, Bd = 0.32, Dd = 0.2
    if (.not.assert(dsys%A, reshape([0.6d0], [1, 1]), tol)) rst = .false.
    if (.not.assert(dsys%B, reshape([0.32d0], [1, 1]), tol)) rst = .false.
    if (.not.assert(dsys%D, reshape([0.2d0], [1, 1]), tol)) rst = .false.
end function

! ------------------------------------------------------------------------------
function test_discretize_round_trip() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    integer :: i
    integer, dimension(2) :: methods
    type(discrete_state_space) :: dsys
    type(state_space) :: back, sys

    ! Initialization
    rst = .true.
    methods = [DISCRETIZE_ZOH, DISCRETIZE_TUSTIN]
    sys = state_space(reshape([-1.0d0, 0.0d0, 2.0d0, -3.0d0], [2, 2]), &
        reshape([1.0d0, 1.0d0], [2, 1]), reshape([1.0d0, 0.0d0], [1, 2]), &
        reshape([0.5d0], [1, 1]))

    ! Test: continuous -> discrete -> continuous recovers the original model
    do i = 1, size(methods)
        dsys = c2d(sys, 0.2d0, methods(i))
        back = d2c(dsys, methods(i))
        if (.not.assert(back%A, sys%A, tol)) rst = .false.
        if (.not.assert(back%B, sys%B, tol)) rst = .false.
        if (.not.assert(back%C, sys%C, tol)) rst = .false.
        if (.not.assert(back%D, sys%D, tol)) rst = .false.
    end do
end function

! ------------------------------------------------------------------------------
end module
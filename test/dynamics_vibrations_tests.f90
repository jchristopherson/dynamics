module dynamics_vibrations_tests
    use iso_fortran_env
    use dynamics
    use fortran_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
function test_q_factor() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: q, ans, zeta

    ! Initialization
    rst = .true.
    call random_number(zeta)
    ans = 1.0d0 / (2.0d0 * zeta)

    ! Test
    q = q_factor(zeta)
    if (.not.assert(q, ans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_q_factor -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_bandwidth() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: df, ans, fn, zeta

    ! Initialization
    rst = .true.
    call random_number(zeta)
    call random_number(fn)
    ans = 2.0d0 * zeta * fn

    ! Test
    df = estimate_bandwidth(fn, zeta)
    if (.not.assert(df, ans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_bandwidth -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_log_decrement() result(rst)
    ! Arguments
    logical :: rst

    ! Variables
    integer(int32), parameter :: n = 2
    real(real64) :: x1, x2, ans, delta

    ! Initialization
    rst = .true.
    call random_number(x1)
    x2 = 0.8d0 * x1
    ans = log(x1 / x2) / n

    ! Test
    delta = logarithmic_decrement(x1, x2, n)
    if (.not.assert(delta, ans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_log_decrement -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_damping_from_decrement() result(rst)
    ! Arguments
    logical :: rst

    ! Variables
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64) :: delta, zeta, ans

    ! Initialization
    rst = .true.
    call random_number(delta)
    ans = delta / sqrt(4.0d0 * pi**2 + delta**2)

    ! Test
    zeta = damping_from_log_decrement(delta)
    if (.not.assert(zeta, ans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_damping_from_decrement -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_find_free_rsp_props() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: fn = 1.0d1
    real(real64), parameter :: zeta = 1.0d-1
    real(real64), parameter :: xo = 0.0d0
    real(real64), parameter :: vo = 1.0d0
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: dt = 1.0d-3
    integer(int32), parameter :: n = 1000

    ! Variables
    integer(int32) :: i
    real(real64) :: t(n), x(n), wn, wd, A, B, fnx, delta, zx, tol

    ! Initialization
    rst = .true.
    t = (/ (i * dt, i = 0, n - 1) /)
    wn = 2.0d0 * pi * fn
    wd = wn * sqrt(1.0d0 - zeta**2)
    A = (vo + zeta * wn * xo) / wd
    B = xo
    x = exp(-zeta * wn * t) * (A * sin(wd * t) + B * cos(wd * t))
    tol = 1.0d-3

    ! Test
    call find_free_response_properties(t, x, delta, fnx)
    zx = damping_from_log_decrement(delta)

    if (.not.assert(zx, zeta, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_find_free_rsp_props -1"
    end if

    tol = abs(wn - wd)
    if (.not.assert(fnx, fn, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_find_free_rsp_props -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_rise_time() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    real(real64) :: wn, zeta, ans, tr

    ! Initialization
    rst = .true.
    call random_number(wn)
    call random_number(zeta)
    ans = (1.0d0 / (wn * sqrt(1.0d0 - zeta**2))) * (pi - &
        atan(sqrt(1.0d0 - zeta**2) / zeta))

    ! Test
    tr = rise_time(wn, zeta)
    if (.not.assert(ans, tr)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_rise_time -1"
    end if
end function

! ------------------------------------------------------------------------------
end module
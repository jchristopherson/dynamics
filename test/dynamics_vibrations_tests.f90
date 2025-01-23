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

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
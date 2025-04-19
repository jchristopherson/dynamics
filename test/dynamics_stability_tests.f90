module dynamics_stability_tests
    use iso_fortran_env
    use dynamics
    use fortran_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
function test_local_stability() result(rst)
    ! Arguments
    logical :: rst

    ! This test utilizes the unforced Duffing's equation of the form:
    ! x" + d x' + a x + b x**3 = 0

    ! Parameters
    real(real64), parameter :: d = 1.0d-1
    real(real64), parameter :: a = 1.0d0
    real(real64), parameter :: b = -5.0d0

    ! Local Variables
    integer(int32) :: flag
    real(real64) :: pt1, pt2, pt3, A1(2, 2), A2(2, 2), A3(2, 2)
    complex(real64) :: ev1(2), ev2(2), ev3(2)

    ! Initialization
    rst = .true.

    ! The equilibrium condition can be determined from the unforced system when
    ! x" = x' = 0.
    !
    ! a x + b x**3 = 0
    !
    ! Solving for x:
    ! x = 0, +/- sqrt(-a / b)
    pt1 = 0.0d0
    pt2 = sqrt(-a / b)
    pt3 = -pt2

    ! The system written as a set of 1st order equations:
    ! F1 = x1 = x'
    ! F2 = x1' = -a x - b x**3 - d x1
    !
    ! The Jacobian matrix:
    !     | dF1 / dx    dF1 / dx1 |   |     0                1 |
    ! A = |                       | = |                        |
    !     | dF2 / dx    dF2 / dx1 |   | -a - 3 b x**2       -d |
    A1 = reshape([0.0d0,  -a - 3.0d0 * b * pt1**2, 1.0d0, -d], [2, 2])
    A2 = reshape([0.0d0,  -a - 3.0d0 * b * pt2**2, 1.0d0, -d], [2, 2])
    A3 = reshape([0.0d0,  -a - 3.0d0 * b * pt3**2, 1.0d0, -d], [2, 2])

    ! Compute the stability condition around pt1
    flag = determine_local_stability(A1, ev = ev1)
    if (flag /= HYPERBOLIC_FIXED_POINT_SINK) then
        rst = .false.
        print "(A)", "TEST: test_local_stability -1"
        print *, ev1
    end if
    
    ! Compute the stability condition around pt2
    flag = determine_local_stability(A2, ev = ev2)
    if (flag /= HYPERBOLIC_FIXED_POINT_SADDLE) then
        rst = .false.
        print "(A)", "TEST: test_local_stability -2"
        print *, ev2
    end if

    ! Compute the stability condition around pt3
    flag = determine_local_stability(A3, ev = ev3)
    if (flag /= HYPERBOLIC_FIXED_POINT_SADDLE) then
        rst = .false.
        print "(A)", "TEST: test_local_stability -3"
        print *, ev3
    end if
end function

! ------------------------------------------------------------------------------
end module
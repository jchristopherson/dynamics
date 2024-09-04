module dynamics_kinematics_tests
    use iso_fortran_env
    use dynamics
    use fortran_test_helper
    implicit none

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

contains
! ------------------------------------------------------------------------------
function test_forward_kinematics() result(rst)
    ! Arguments
    logical :: rst

    ! Mechanism: 3R PUMA

    ! Link Lengths
    real(real64), parameter :: L1 = 1.25d0
    real(real64), parameter :: L2 = 5.5d0
    real(real64), parameter :: L3 = 2.0d0

    ! DH Parameters
    real(real64), parameter :: a1 = 0.d0
    real(real64), parameter :: a2 = L2;
    real(real64), parameter :: a3 = L3;

    real(real64), parameter :: alpha1 = 0.5d0 * pi
    real(real64), parameter :: alpha2 = 0.0d0
    real(real64), parameter :: alpha3 = 0.0d0
    
    real(real64), parameter :: theta1 = 0.0d0
    real(real64), parameter :: theta2 = 0.5d0 * pi
    real(real64), parameter :: theta3 = 0.25d0 * pi

    real(real64), parameter :: d1 = 0.0d0
    real(real64), parameter :: d2 = -L1
    real(real64), parameter :: d3 = 0.0d0

    ! Local Variables
    real(real64), dimension(4,4) :: T1, T2, T3, Tref, T

    ! Initialization
    rst = .true.

    ! Compute individual transformation matrices for each link
    T1 = dh_matrix(alpha1, a1, theta1, d1)
    T2 = dh_matrix(alpha2, a2, theta2, d2)
    T3 = dh_matrix(alpha3, a3, theta3, d3)

    ! Compute the reference matrix
    Tref = matmul(T1, matmul(T2, T3))

    ! Now, compute using the forward kinematics routine
    T = dh_forward_kinematics(T1, T2, T3)

    ! Test
    if (.not.assert(Tref, T)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_forward_kinematics -1"
    end if
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module
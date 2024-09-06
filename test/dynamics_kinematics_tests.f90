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

    ! Assemble via arrays
    T = dh_forward_kinematics( &
        [alpha1, alpha2, alpha3], &
        [a1, a2, a3], &
        [theta1, theta2, theta3], &
        [d1, d2, d3] &
    )
    if (.not.assert(Tref, T)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_forward_kinematics -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_inverse_kinematics() result(rst)
    ! Arguments
    logical :: rst

    ! Tolerance
    real(real64), parameter :: tol = 1.0d-4

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
    real(real64) :: T(4, 4), ans(3), q(3), qo(3), Tc(4, 4), constraints(6)
    procedure(vecfcn), pointer :: mdl

    ! Initialization
    rst = .true.
    ans = [theta1, theta2, theta3]

    ! Solve the forward kinematics problem to define the constraints
    T = dh_forward_kinematics( &
        [alpha1, alpha2, alpha3], &
        [a1, a2, a3], &
        [theta1, theta2, theta3], &
        [d1, d2, d3] &
    )
    constraints = [ &
        T(1,4), T(2,4), T(3, 4), &
        T(1,1), T(2,1), T(3,1) &
    ]

    ! Define an initial guess
    qo = 0.0d0

    ! Solve the inverse problem
    mdl => inverse_test_kinematics_equation
    q = solve_inverse_kinematics(mdl, qo, constraints)

    ! Verify the transformation matrix
    Tc = dh_forward_kinematics( &
        [alpha1, alpha2, alpha3], &
        [a1, a2, a3], &
        q, &
        [d1, d2, d3] &
    )

    ! Test
    if (.not.assert(q, ans, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_inverse_kinematics -1"
    end if
end function

! -----
subroutine inverse_test_kinematics_equation(q, f)
    ! Arguments
    real(real64), intent(in), dimension(:) :: q
    real(real64), intent(out), dimension(:) :: f

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

    real(real64), parameter :: d1 = 0.0d0
    real(real64), parameter :: d2 = -L1
    real(real64), parameter :: d3 = 0.0d0

    ! Local Variables
    real(real64) :: T(4, 4)

    ! Define the location of the end-effector
    T = dh_forward_kinematics( &
        [alpha1, alpha2, alpha3], &
        [a1, a2, a3], &
        q, &
        [d1, d2, d3] &
    )

    ! Define the position of the end effector
    f(1:3) = T(1:3,4)

    ! Define its orientation
    f(4:6) = T(1:3,1)
end subroutine

! ------------------------------------------------------------------------------
end module
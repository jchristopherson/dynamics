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
subroutine inverse_test_kinematics_equation(q, f, args)
    ! Arguments
    real(real64), intent(in), dimension(:) :: q
    real(real64), intent(out), dimension(:) :: f
    class(*), intent(inout), optional :: args

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
function test_jacobian() result(rst)
    ! This is Example 214 from Jazar's kinematics text.

    ! Arguments
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-8
    real(real64), parameter :: zi_1(4) = [0.0d0, 0.0d0, 1.0d0, 0.0d0]

    ! Local Variables
    real(real64) :: l0, theta1, theta2, d3, alpha(3), theta(3), a(3), d(3), &
        c1(6), c2(6), c3(6), Jans(6, 3), T1(4, 4), T2(4, 4), T3(4, 4), &
        di_1(3), ki_1(4), T1ans(4, 4), T2ans(4, 4), T3ans(4, 4), T2a(4, 4), &
        T3a(4, 4), T(4, 4), k1ans(3), d1ans(3), k2ans(3), R(3, 3), J(6, 3)
    integer(int32) :: jtypes(3)

    ! Initialization
    rst = .true.
    call random_number(l0)
    call random_number(theta1)
    call random_number(theta2)
    call random_number(d3)
    
    alpha = [-0.5d0 * pi, 0.5d0 * pi, 0.0d0]
    theta = [theta1, theta2, 0.0d0]
    a = [0.0d0, 0.0d0, 0.0d0]
    d = [l0, 0.0d0, d3]

    jtypes = [ &
        REVOLUTE_JOINT, &
        REVOLUTE_JOINT, &
        PRISMATIC_JOINT &
    ]

    ! Compute the transformation matrices
    T1 = dh_matrix(alpha(1), a(1), theta(1), d(1))
    T2a = dh_matrix(alpha(2), a(2), theta(2), d(2))
    T2 = matmul(T1, T2a)
    T3a = dh_matrix(alpha(3), a(3), theta(3), d(3))
    T3 = matmul(T2, T3a)

    ! The Jacobian generating vectors
    Jans(:,1) = [ &
        -d3 * sin(theta1) * sin(theta2), &
        d3 * cos(theta1) * sin(theta2), &
        0.0d0, &
        0.0d0, &
        0.0d0, &
        1.0d0 &
    ]
    Jans(:,2) = [ &
        d3 * cos(theta1) * cos(theta2), &
        d3 * cos(theta2) * sin(theta1), &
        -d3 * sin(theta2), &
        -sin(theta1), &
        cos(theta1), &
        0.0d0 &
    ]
    Jans(:,3) = [ &
        cos(theta1) * sin(theta2), &
        sin(theta1) * sin(theta2), &
        cos(theta2), &
        0.0d0, &
        0.0d0, &
        0.0d0 &
    ]

    ! Verify the transformation matrices are correct
    T1ans = reshape([ &
        cos(theta1), sin(theta1), 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, -1.0d0, 0.0d0, &
        -sin(theta1), cos(theta1), 0.0d0, 0.0d0, &
        0.0d0, 0.0d0, l0, 1.0d0 &
    ], [4, 4])
    T2ans = reshape([ &
        cos(theta1) * cos(theta2), cos(theta2) * sin(theta1), -sin(theta2), 0.0d0, &
        -sin(theta1), cos(theta1), 0.0d0, 0.0d0, &
        cos(theta1) * sin(theta2), sin(theta1) * sin(theta2), cos(theta2), 0.0d0, &
        0.0d0, 0.0d0, l0, 1.0d0 &
    ], [4, 4])
    T3ans = matmul(T2ans, dh_translate_z(d3))

    if (.not.assert(T1, T1ans, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_jacobian -1"
    end if
    if (.not.assert(T2, T2ans, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_jacobian -2"
    end if
    if (.not.assert(T3, T3ans, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_jacobian -3"
    end if

    ! Compute each Jacobian generating vector
    di_1 = T3(1:3,4)
    d1ans = [&
        d3 * cos(theta1) * sin(theta2), &
        d3 * sin(theta1) * sin(theta2), &
        l0 + d3 * cos(theta2) &
    ]
    if (.not.assert(di_1, d1ans, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_jacobian -4"
    end if
    ki_1 = zi_1
    k1ans = [0.0d0, 0.0d0, 1.0d0]
    if (.not.assert(ki_1(1:3), k1ans, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_jacobian -5"
    end if
    R = reshape( &
        [1.0d0, 0.0d0, 0.0d0, &
        0.0d0, 1.0d0, 0.0d0, &
        0.0d0, 0.0d0, 1.0d0], & 
        [3, 3])
    c1 = jacobian_generating_vector(di_1, ki_1(1:3), R, jtypes(1))
    if (.not.assert(c1, Jans(:,1), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_jacobian -6"
    end if

    ! -----
    T = matmul(T2a, T3a)
    di_1 = T(1:3,4)
    ki_1 = matmul(T1, zi_1)
    k2ans = [ &
        -sin(theta1), &
        cos(theta1), &
        0.0d0 &
    ]
    if (.not.assert(ki_1(1:3), k2ans, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_jacobian -7"
    end if
    R = T1(1:3,1:3)
    c2 = jacobian_generating_vector(di_1, ki_1(1:3), R, jtypes(2))
    if (.not.assert(c2, Jans(:,2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_jacobian -8"
    end if

    ! -----
    T = T3a
    di_1 = T(1:3,4)
    ki_1 = matmul(T2, zi_1)
    R = T2(1:3,1:3)
    c3 = jacobian_generating_vector(di_1, ki_1(1:3), R, jtypes(3))
    if (.not.assert(c3, Jans(:,3), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_jacobian -9"
    end if

    ! Test the entire Jacobian
    J = dh_jacobian(alpha, a, theta, d, jtypes)
    if (.not.assert(J, Jans, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_jacobian -10"
    end if
end function

! ------------------------------------------------------------------------------
function test_velocity_matrix() result(rst)
    ! Arguments
    logical :: rst

    ! Variables
    real(real64) :: omega(3), v(3), x(3), T(4, 4), ans(4, 4), w(3, 3)

    ! Initialization
    rst = .true.
    call random_number(omega)
    call random_number(v)
    call random_number(x)
    w = to_skew_symmetric(omega)

    ! Compute the solution
    ans(1:3,1:3) = w
    ans(1:3,4) = v - matmul(w, x)
    ans(4,:) = 0.0d0

    ! Test
    T = velocity_transform(omega, v, x)
    if (.not.assert(ans, T)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_velocity_matrix -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_acceleration_matrix() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: alpha(3), a(3), omega(3), x(3), T(4, 4), al(3,3), w(3,3), &
        ans(4, 4)

    ! Initialization
    rst = .true.
    call random_number(alpha)
    call random_number(a)
    call random_number(omega)
    call random_number(x)
    al = to_skew_symmetric(alpha)
    w = to_skew_symmetric(omega)

    ! Compute the answer
    ans(1:3,1:3) = al - matmul(w, transpose(w))
    ans(1:3,4) = a - matmul(ans(1:3,1:3), x)
    ans(4,:) = 0.0d0

    ! Test
    T = acceleration_transform(alpha, omega, a, x)
    if (.not.assert(T, ans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_acceleration_matrix -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_skew_symmetric() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: x(3), xt(3,3), ans(3,3)

    ! Initialization
    rst = .true.
    call random_number(x)
    ans = reshape([ &
        0.0d0, x(3), -x(2), &
        -x(3), 0.0d0, x(1), &
        x(2), -x(1), 0.0d0 &
    ], [3, 3])

    ! Test
    xt = to_skew_symmetric(x)
    if (.not.assert(xt, ans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_skew_symmetric -1"
    end if
end function

! ------------------------------------------------------------------------------
end module
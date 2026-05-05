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
! Parallel joint axes
function test_define_link_csys_1() result(rst)
    logical :: rst
    real(real64) :: xim1(3), zim1(3), zi(3), rim1(3), ri(3), u(3), &
        ians(3), jans(3), kans(3)
    type(coordinate_system) :: csys

    ! Initialization
    rst = .true.
    xim1 = [1.0d0, 0.0d0, 0.0d0]
    zim1 = [0.0d0, 0.0d0, 1.0d0]
    zi = [0.0d0, 0.0d0, 1.0d0]
    rim1 = [0.0d0, 0.0d0, 0.0d0]
    ri = [1.0d0, 1.0d0, 0.0d0]
    u = ri - rim1
    ians = u / norm2(u)
    kans = zi
    jans = cross_product(kans, ians)

    ! Build the coordinate system
    csys = coordinate_system(xim1, zim1, zi, rim1, ri)

    ! Tests
    if (.not.assert(csys%i, ians)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_1 -1"
    end if

    if (.not.assert(csys%j, jans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_1 -2"
    end if

    if (.not.assert(csys%k, kans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_1 -3"
    end if

    if (.not.assert(csys%origin, u)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_1 -4"
    end if
end function

! ------------------------------------------------------------------------------
! Skew, non-intersecting joint axes
function test_define_link_csys_2() result(rst)
    logical :: rst
    real(real64) :: xim1(3), zim1(3), zi(3), rim1(3), ri(3), u(3), &
        ians(3), jans(3), kans(3)
    type(coordinate_system) :: csys

    ! Initialization
    rst = .true.
    xim1 = [1.0d0, 0.0d0, 0.0d0]
    zim1 = [0.0d0, 0.0d0, 1.0d0]
    zi = [1.0d0, 0.0d0, 0.0d0]
    rim1 = [0.0d0, 0.0d0, 0.0d0]
    ri = [1.0d0, 1.0d0, 0.0d0]
    u = [0.0d0, 1.0d0, 0.0d0]
    ians = [0.0d0, 1.0d0, 0.0d0]
    kans = zi
    jans = cross_product(kans, ians)

    ! Build the coordinate system
    csys = coordinate_system(xim1, zim1, zi, rim1, ri)

    ! Tests
    if (.not.assert(csys%i, ians)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_2 -1"
    end if

    if (.not.assert(csys%j, jans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_2 -2"
    end if

    if (.not.assert(csys%k, kans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_2 -3"
    end if

    if (.not.assert(csys%origin, u)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_2 -4"
    end if
end function

! ------------------------------------------------------------------------------
! Intersecting joint axes
function test_define_link_csys_3() result(rst)
    logical :: rst
    real(real64) :: xim1(3), zim1(3), zi(3), rim1(3), ri(3), u(3), &
        ians(3), jans(3), kans(3)
    type(coordinate_system) :: csys

    ! Initialization
    rst = .true.
    xim1 = [1.0d0, 0.0d0, 0.0d0]
    zim1 = [0.0d0, 0.0d0, 1.0d0]
    zi = [0.0d0, 1.0d0, 0.0d0]
    rim1 = [1.0d0, 1.0d0, 1.0d0]
    ri = [1.0d0, 1.0d0, -1.0d0]
    u = ri
    ians = cross_product(zim1, zi)
    kans = zi
    jans = cross_product(kans, ians)

    ! Build the coordinate system
    csys = coordinate_system(xim1, zim1, zi, rim1, ri)

    ! Tests
    if (.not.assert(csys%i, ians)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_3 -1"
    end if

    if (.not.assert(csys%j, jans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_3 -2"
    end if

    if (.not.assert(csys%k, kans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_3 -3"
    end if

    if (.not.assert(csys%origin, u)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_3 -4"
    end if
end function

! ------------------------------------------------------------------------------
! Colinear joint axes
function test_define_link_csys_4() result(rst)
    logical :: rst
    real(real64) :: xim1(3), zim1(3), zi(3), rim1(3), ri(3), u(3), &
        ians(3), jans(3), kans(3)
    type(coordinate_system) :: csys

    ! Initialization
    rst = .true.
    xim1 = [1.0d0, 0.0d0, 0.0d0]
    zim1 = [0.0d0, 0.0d0, 1.0d0]
    zi = [0.0d0, 0.0d0, 1.0d0]
    rim1 = [0.0d0, 0.0d0, 0.0d0]
    ri = [0.0d0, 0.0d0, 1.0d0]
    u = ri
    ians = xim1
    kans = zi
    jans = cross_product(kans, ians)

    ! Build the coordinate system
    csys = coordinate_system(xim1, zim1, zi, rim1, ri)

    ! Tests
    if (.not.assert(csys%i, ians)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_4 -1"
    end if

    if (.not.assert(csys%j, jans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_4 -2"
    end if

    if (.not.assert(csys%k, kans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_4 -3"
    end if

    if (.not.assert(csys%origin, u)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_define_link_csys_4 -4"
    end if
end function

! ------------------------------------------------------------------------------
! 2D Serial Manipulator
!
! Jazar's Text
! Example 126
! Page 204
function test_dh_table_1() result(rst)
    logical :: rst
    real(real64) :: l1, l2, l3, t1, t2, t3, R1(3,3), R2(3,3), R3(3,3), R(3,3), &
        i(3), j(3), p1(3), p2(3), p3(3)
    type(coordinate_system) :: csys(4)
    type(dh_table) :: tbl

    ! Parameters
    real(real64), parameter :: ib(3) = [1.0d0, 0.0d0, 0.0d0]
    real(real64), parameter :: jb(3) = [0.0d0, 1.0d0, 0.0d0]
    real(real64), parameter :: k(3) = [0.0d0, 0.0d0, 1.0d0]

    ! Initialization
    rst = .true.
    call random_number(l1)
    call random_number(l2)
    call random_number(l3)
    call random_number(t1)
    call random_number(t2)
    call random_number(t3)
    t2 = -t2    ! flip sign on an angle as a sanity check

    ! Define the coordinate systems
    R1 = rotate_z(t1)
    R2 = rotate_z(t2)
    R3 = rotate_z(t3)

    ! CSYS 0
    csys(1) = coordinate_system(ib, jb, k, [0.0d0, 0.0d0, 0.0d0])

    ! CSYS 1
    i = matmul(R1, ib)
    j = matmul(R1, jb)
    p1 = matmul(R1, l1 * ib)
    csys(2) = coordinate_system(i, j, k, p1)

    ! CSYS 2
    R = matmul(R1, R2)
    i = matmul(R, ib)
    j = matmul(R, jb)
    p2 = matmul(R, l2 * ib) + p1
    csys(3) = coordinate_system(i, j, k, p2)

    ! CSYS 3
    R = matmul(R, R3)
    i = matmul(R, ib)
    j = matmul(R, jb)
    p3 = matmul(R, l3 * ib) + p2
    csys(4) = coordinate_system(i, j, k, p3)

    ! Build the table
    tbl = dh_table(csys)

    ! Tests
    if (size(tbl%parameters) /= 3) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1"
        return
    end if

    if (.not.assert(tbl%parameters(1)%link_length, l1)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -1a"
    end if
    if (.not.assert(tbl%parameters(1)%link_twist, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -1b"
    end if
    if (.not.assert(tbl%parameters(1)%link_offset, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -1c"
    end if
    if (.not.assert(tbl%parameters(1)%joint_angle, t1)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -1d"
    end if

    if (.not.assert(tbl%parameters(2)%link_length, l2)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -2a"
    end if
    if (.not.assert(tbl%parameters(2)%link_twist, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -2b"
    end if
    if (.not.assert(tbl%parameters(2)%link_offset, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -2c"
    end if
    if (.not.assert(tbl%parameters(2)%joint_angle, t2)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -2d"
    end if

    if (.not.assert(tbl%parameters(3)%link_length, l3)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -3a"
    end if
    if (.not.assert(tbl%parameters(3)%link_twist, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -3b"
    end if
    if (.not.assert(tbl%parameters(3)%link_offset, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -3c"
    end if
    if (.not.assert(tbl%parameters(3)%joint_angle, t3)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -3d"
    end if
end function

! ------------------------------------------------------------------------------
! 3D PUMA Manipulator
!
! Jazar's Text
! Example 127
! Page 204, 205
function test_dh_table_2() result(rst)
    logical :: rst
    real(real64) :: l1, l2, l3, p1(3), p2(3), p3(3)
    type(coordinate_system) :: csys(4)
    type(dh_table) :: tbl

    ! Parameters
    real(real64), parameter :: ib(3) = [1.0d0, 0.0d0, 0.0d0]
    real(real64), parameter :: jb(3) = [0.0d0, 1.0d0, 0.0d0]
    real(real64), parameter :: kb(3) = [0.0d0, 0.0d0, 1.0d0]
    real(real64), parameter :: half_pi = acos(0.0d0)

    ! Initialization
    rst = .true.
    ! call random_number(l1)
    ! call random_number(l2)
    ! call random_number(l3)
    l1 = 0.5d0
    l2 = 1.5d0
    l3 = 1.0d0
    p1 = [0.0d0, 0.0d0, 0.0d0]
    p2 = [l2, l1, 0.0d0]
    p3 = [0.0d0, 0.0d0, l3] + p2

    ! CSYS 0
    csys(1) = coordinate_system(ib, jb, kb, [0.0d0, 0.0d0, 0.0d0])

    ! CSYS 1
    csys(2) = coordinate_system(ib, kb, -jb, p1)

    ! CSYS 2
    csys(3) = coordinate_system(ib, kb, -jb, p2)

    ! CSYS 3
    csys(4) = coordinate_system(kb, -ib, -jb, p3)

    ! Build the table
    tbl = dh_table(csys)

    ! Test
    if (size(tbl%parameters) /= 3) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1"
        return
    end if

    if (.not.assert(tbl%parameters(1)%link_length, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -1a"
    end if
    if (.not.assert(tbl%parameters(1)%link_twist, half_pi)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -1b"
    end if
    if (.not.assert(tbl%parameters(1)%link_offset, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -1c"
    end if
    if (.not.assert(tbl%parameters(1)%joint_angle, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_1 -1d"
    end if

    if (.not.assert(tbl%parameters(2)%link_length, l2)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_2 -2a"
    end if
    if (.not.assert(tbl%parameters(2)%link_twist, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_2 -2b"
    end if
    if (.not.assert(tbl%parameters(2)%link_offset, l1)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_2 -2c"
    end if
    if (.not.assert(tbl%parameters(2)%joint_angle, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_2 -2d"
    end if

    if (.not.assert(tbl%parameters(3)%link_length, l3)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_2 -3a"
    end if
    if (.not.assert(tbl%parameters(3)%link_twist, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_2 -3b"
    end if
    if (.not.assert(tbl%parameters(3)%link_offset, 0.0d0)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_2 -3c"
    end if
    if (.not.assert(tbl%parameters(3)%joint_angle, half_pi)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_2 -3d"
    end if
end function

! ------------------------------------------------------------------------------
! Stanford Arm 
!
! Jazar's Text
! Example 128
! Page 205, 206
function test_dh_table_3() result(rst)
    logical :: rst
    real(real64) :: l1, l2, l3, t1, t2, t4, t5, t6, l6
    real(real64), dimension(4) :: p1, p2, p3, p4, p5, p6
    real(real64), dimension(4) :: i1, j1, k1
    real(real64), dimension(4) :: i2, j2, k2
    real(real64), dimension(4) :: i3, j3, k3
    real(real64), dimension(4) :: i4, j4, k4
    real(real64), dimension(4) :: i5, j5, k5
    real(real64), dimension(4) :: i6, j6, k6
    real(real64) :: a1, alpha1, d1, theta1
    real(real64) :: a2, alpha2, d2, theta2
    real(real64) :: a3, alpha3, d3, theta3
    real(real64) :: a4, alpha4, d4, theta4
    real(real64) :: a5, alpha5, d5, theta5
    real(real64) :: a6, alpha6, d6, theta6
    real(real64), dimension(4,4) :: T01, T12, T23, T34, T45, T56, T, Tt
    real(real64) :: R1(3,3)
    type(coordinate_system) :: csys(7)
    type(dh_table) :: tbl

    ! Parameters
    real(real64), parameter :: ib(4) = [1.0d0, 0.0d0, 0.0d0, 0.0d0]
    real(real64), parameter :: jb(4) = [0.0d0, 1.0d0, 0.0d0, 0.0d0]
    real(real64), parameter :: kb(4) = [0.0d0, 0.0d0, 1.0d0, 0.0d0]
    real(real64), parameter :: origin(4) = [0.0d0, 0.0d0, 0.0d0, 1.0d0]
    real(real64), parameter :: half_pi = acos(0.0d0)

    ! Initialization
    rst = .true.
    l1 = 1.5d0
    l2 = 0.25d0
    l3 = 5.0d0
    l6 = 2.5d0
    t1 = 0.25d0
    t2 = -0.1d0
    t4 = 0.5d0
    t5 = -0.3d0
    t6 = 0.05d0

    ! Answers
    a1 = 0.0d0
    alpha1 = -half_pi
    d1 = l1
    theta1 = t1

    a2 = 0.0d0
    alpha2 = half_pi
    d2 = l2
    theta2 = t2

    a3 = 0.0d0
    alpha3 = 0.0d0
    d3 = l3
    theta3 = 0.0d0

    a4 = 0.0d0
    alpha4 = -half_pi
    d4 = 0.0d0
    theta4 = t4

    a5 = 0.0d0
    alpha5 = half_pi
    d5 = 0.0d0
    theta5 = t5

    a6 = 0.0d0
    alpha6 = 0.0d0
    d6 = l6
    theta6 = t6

    ! Joint Positions & Orientations
    T01 = dh_matrix(alpha1, a1, theta1, d1)
    T12 = dh_matrix(alpha2, a2, theta2, d2)
    T23 = dh_matrix(alpha3, a3, theta3, d3)
    T34 = dh_matrix(alpha4, a4, theta4, d4)
    T45 = dh_matrix(alpha5, a5, theta5, d5)
    T56 = dh_matrix(alpha6, a6, theta6, d6)

    p1 = matmul(T01, origin)
    i1 = matmul(T01, ib)
    j1 = matmul(T01, jb)
    k1 = matmul(T01, kb)

    T = matmul(T01, T12)
    p2 = matmul(T, origin)
    i2 = matmul(T, ib)
    j2 = matmul(T, jb)
    k2 = matmul(T, kb)

    T = matmul(T, T23)
    p3 = matmul(T, origin)
    i3 = matmul(T, ib)
    j3 = matmul(T, jb)
    k3 = matmul(T, kb)

    T = matmul(T, T34)
    p4 = matmul(T, origin)
    i4 = matmul(T, ib)
    j4 = matmul(T, jb)
    k4 = matmul(T, kb)

    T = matmul(T, T45)
    p5 = matmul(T, origin)
    i5 = matmul(T, ib)
    j5 = matmul(T, jb)
    k5 = matmul(T, kb)

    T = matmul(T, T56)
    p6 = matmul(T, origin)
    i6 = matmul(T, ib)
    j6 = matmul(T, jb)
    k6 = matmul(T, kb)

    ! Coordinate Systems
    csys(1) = coordinate_system(ib(1:3), jb(1:3), kb(1:3), origin(1:3))
    csys(2) = coordinate_system(i1(1:3), j1(1:3), k1(1:3), p1(1:3))
    csys(3) = coordinate_system(i2(1:3), j2(1:3), k2(1:3), p2(1:3))
    csys(4) = coordinate_system(i3(1:3), j3(1:3), k3(1:3), p3(1:3))
    csys(5) = coordinate_system(i4(1:3), j4(1:3), k4(1:3), p4(1:3))
    csys(6) = coordinate_system(i5(1:3), j5(1:3), k5(1:3), p5(1:3))
    csys(7) = coordinate_system(i6(1:3), j6(1:3), k6(1:3), p6(1:3))

    ! Build the table
    tbl = dh_table(csys)

    ! Tests
    if (size(tbl%parameters) /= 6) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3"
        return
    end if

    if (.not.assert(tbl%parameters(1)%link_length, a1)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -1a"
    end if
    if (.not.assert(tbl%parameters(1)%link_twist, alpha1)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -1b"
    end if
    if (.not.assert(tbl%parameters(1)%link_offset, d1)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -1c"
    end if
    if (.not.assert(tbl%parameters(1)%joint_angle, theta1)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -1d"
    end if

    if (.not.assert(tbl%parameters(2)%link_length, a2)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -2a"
    end if
    if (.not.assert(tbl%parameters(2)%link_twist, alpha2)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -2b"
    end if
    if (.not.assert(tbl%parameters(2)%link_offset, d2)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -2c"
    end if
    if (.not.assert(tbl%parameters(2)%joint_angle, theta2)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -2d"
    end if

    if (.not.assert(tbl%parameters(3)%link_length, a3)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -3a"
    end if
    if (.not.assert(tbl%parameters(3)%link_twist, alpha3)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -3b"
    end if
    if (.not.assert(tbl%parameters(3)%link_offset, d3)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -3c"
    end if
    if (.not.assert(tbl%parameters(3)%joint_angle, theta3)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -3d"
    end if

    if (.not.assert(tbl%parameters(4)%link_length, a4)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -4a"
    end if
    if (.not.assert(tbl%parameters(4)%link_twist, alpha4)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -4b"
    end if
    if (.not.assert(tbl%parameters(4)%link_offset, d4)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -4c"
    end if
    if (.not.assert(tbl%parameters(4)%joint_angle, theta4)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -4d"
    end if

    if (.not.assert(tbl%parameters(5)%link_length, a5)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -5a"
    end if
    if (.not.assert(tbl%parameters(5)%link_twist, alpha5)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -5b"
    end if
    if (.not.assert(tbl%parameters(5)%link_offset, d5)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -5c"
    end if
    if (.not.assert(tbl%parameters(5)%joint_angle, theta5)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -5d"
    end if

    if (.not.assert(tbl%parameters(6)%link_length, a6)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -6a"
    end if
    if (.not.assert(tbl%parameters(6)%link_twist, alpha6)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -6b"
    end if
    if (.not.assert(tbl%parameters(6)%link_offset, d6)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -6c"
    end if
    if (.not.assert(tbl%parameters(6)%joint_angle, theta6)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -6d"
    end if

    ! Test constructing the matrix from the table
    Tt = dh_forward_kinematics(tbl)
    if (.not.assert(T, Tt)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_dh_table_3 -7"
    end if
end function

! ------------------------------------------------------------------------------
end module
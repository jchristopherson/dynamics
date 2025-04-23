module dynamics_transfer_function_tests
    use iso_fortran_env
    use dynamics
    use fortran_test_helper
    implicit none

contains
! ------------------------------------------------------------------------------
function test_tf_evaluate() result(rst)
    logical :: rst

    ! Parameters
    complex(real64), parameter :: j = (0.0d0, 1.0d0)
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    complex(real64) :: s, ans, h1, h2
    real(real64) :: omega, m, b, k
    type(transfer_function) :: tf

    ! Initialization
    rst = .true.
    call random_number(omega)
    omega = 1.0d1 * omega
    call random_number(m)
    call random_number(b)
    call random_number(k)
    s = j * omega

    ! Set up the transfer function
    call tf%initialize([1.0d0], [k, b, m])

    ! Compute the solution
    ans = 1.0d0 / (m * s**2 + b * s + k)

    ! Compute using s
    h1 = tf%evaluate(s)
    if (.not.assert(h1, ans, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_tf_evaluate -1"
    end if

    ! Compute using omega
    h2 = tf%evaluate(omega)
    if (.not.assert(h2, ans, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_tf_evaluate -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_tf_multiply() result(rst)
    logical :: rst

    ! Local Variables
    real(real64) :: rc, m, b, k, ans1(1), ans2(4)
    real(real64), allocatable, dimension(:) :: numer, denom
    type(transfer_function) :: tf, tf1, tf2

    ! Initialization
    rst = .true.
    call random_number(rc)
    call random_number(m)
    call random_number(b)
    call random_number(k)
    
    ! Set up the transfer functions
    call tf1%initialize([1.0d0], [k, b, m])
    call tf2%initialize([1.0d0], [1.0d0, rc])

    ! Define the answer
    ans1 = [1.0d0]
    ans2 = [k, rc * k + b, m + rc * b, rc * m]

    ! Compute the solution
    tf = tf1 * tf2

    ! Test
    numer = tf%Y%get_all()
    denom = tf%X%get_all()
    if (.not.assert(numer, ans1)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_tf_multiply -1"
    end if
    if (.not.assert(denom, ans2)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_tf_multiply -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_tf_poles_zeros() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-8

    ! Local Variables
    complex(real64), allocatable, dimension(:) :: pAns, zAns, poles, zeros
    real(real64) :: m, b, k
    type(transfer_function) :: tf

    ! Initialization
    rst = .true.
    call random_number(m)
    call random_number(b)
    call random_number(k)
    call tf%initialize([k, b], [k, b, m])

    ! Compute the solution
    pAns = tf%X%roots()
    zAns = tf%Y%roots()

    ! Poles Test
    poles = tf%poles()
    if (.not.assert(poles, pAns, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_tf_poles_zeros -1"
    end if

    ! Zeros Test
    zeros = tf%zeros()
    if (.not.assert(zeros, zAns, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_tf_poles_zeros -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_ccf_form_conversion() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-8

    ! Variables
    integer(int32), parameter :: n = 4
    integer(int32) :: i
    real(real64) :: A(n,n), B(n,1), C(1,n), D(1,1), x(n + 1), y(n)
    type(transfer_function) :: tf
    type(state_space) :: ss

    ! Initialization
    rst = .true.
    call random_number(x)
    call random_number(y)
    A = 0.0d0
    B = 0.0d0
    C = 0.0d0
    D = 0.0d0
    call tf%initialize(y, x)

    ! Define the solution
    y = y / x(n + 1)
    x = x / x(n + 1)
    A(n,:) = -x(1:n)
    B(n,1) = 1.0d0
    do i = 1, n - 1
        A(i,i+1) = 1.0d0
    end do
    C(1,:) = y

    ! Compute the solution
    ss = tf%to_ccf_state_space()

    ! Tests
    if (.not.assert(A, ss%A, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_ccf_form_conversion -1"
    end if

    if (.not.assert(B, ss%B, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_ccf_form_conversion -2"
    end if

    if (.not.assert(C, ss%C, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_ccf_form_conversion -3"
    end if

    if (.not.assert(D, ss%D, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_ccf_form_conversion -4"
    end if
end function

! ------------------------------------------------------------------------------
function test_ocf_form_conversion() result(rst)
    logical :: rst

    ! Parameters
    real(real64), parameter :: tol = 1.0d-8

    ! Variables
    integer(int32), parameter :: n = 4
    integer(int32) :: i
    real(real64) :: A(n,n), B(n,1), C(1,n), D(1,1), x(n + 1), y(n)
    type(transfer_function) :: tf
    type(state_space) :: ss

    ! Initialization
    rst = .true.
    call random_number(x)
    call random_number(y)
    A = 0.0d0
    B = 0.0d0
    C = 0.0d0
    D = 0.0d0
    call tf%initialize(y, x)

    ! Define the solution
    y = y / x(n + 1)
    x = x / x(n + 1)
    A(:,n) = -x(1:n)
    do i = 1, n - 1
        A(i+1,i) = 1.0d0
    end do
    B(:,1) = y
    C(1,n) = 1.0d0

    ! Compute the solution
    ss = tf%to_ocf_state_space()

    ! Tests
    if (.not.assert(A, ss%A, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_ocf_form_conversion -1"
    end if

    if (.not.assert(B, ss%B, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_ocf_form_conversion -2"
    end if

    if (.not.assert(C, ss%C, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_ocf_form_conversion -3"
    end if

    if (.not.assert(D, ss%D, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_ocf_form_conversion -4"
    end if
end function

! ------------------------------------------------------------------------------
end module
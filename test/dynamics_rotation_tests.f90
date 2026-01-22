module dynamics_rotation_tests
    use iso_fortran_env
    use fortran_test_helper
    use dynamics_rotation
    implicit none

contains
! ------------------------------------------------------------------------------
function test_quaternion_init_1() result(rst)
    ! Variables
    logical :: rst
    type(quaternion) :: q
    real(real64) :: x(4)

    ! Initialization
    rst = .true.
    call random_number(x)
    q = quaternion(x)

    ! Test
    if (.not.assert(x(1), q%w)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_1 - 1"
    end if

    if (.not.assert(x(2), q%x)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_1 - 2"
    end if

    if (.not.assert(x(3), q%y)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_1 - 3"
    end if

    if (.not.assert(x(4), q%z)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_1 - 4"
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_init_2() result(rst)
    ! Variables
    logical :: rst
    type(quaternion) :: q
    real(real64) :: angle, axis(3), w, x, y, z

    ! Initialization
    rst = .true.
    call random_number(angle)
    call random_number(axis)
    w = cos(0.5d0 * angle)
    x = sin(0.5d0 * angle) * axis(1)
    y = sin(0.5d0 * angle) * axis(2)
    z = sin(0.5d0 * angle) * axis(3)
    q = quaternion(angle, axis)

    ! Test
    if (.not.assert(w, q%w)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_2 - 1"
    end if

    if (.not.assert(x, q%x)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_2 - 2"
    end if

    if (.not.assert(y, q%y)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_2 - 3"
    end if

    if (.not.assert(z, q%z)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_2 - 4"
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_init_3() result(rst)
    ! Variables
    logical :: rst
    type(quaternion) :: qx, qy, qz, px, py, pz
    real(real64) :: ax, ay, az, axisX(3), axisY(3), axisZ(3), &
        Rx(3,3), Ry(3,3), Rz(3,3)

    ! Initialization
    rst = .true.
    call random_number(ax)
    call random_number(ay)
    call random_number(az)
    axisX = [1.0d0, 0.0d0, 0.0d0]
    axisY = [0.0d0, 1.0d0, 0.0d0]
    axisZ = [0.0d0, 0.0d0, 1.0d0]
    Rx = rotate_x(ax)
    Ry = rotate_y(ay)
    Rz = rotate_z(az)
    qx = quaternion(ax, axisX)
    qy = quaternion(ay, axisY)
    qz = quaternion(az, axisZ)
    px = quaternion(Rx)
    py = quaternion(Ry)
    pz = quaternion(Rz)

    ! Normalize the quaternions
    call qx%normalize()
    call qy%normalize()
    call qz%normalize()
    call px%normalize()
    call py%normalize()
    call pz%normalize()

    ! Test
    if (.not.assert(qx%w, px%w)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_3 - 1"
    end if

    if (.not.assert(qx%x, px%x)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_3 - 2"
    end if

    if (.not.assert(qx%y, px%y)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_3 - 3"
    end if

    if (.not.assert(qx%z, px%z)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_3 - 4"
    end if

    ! -----
    if (.not.assert(qy%w, py%w)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_3 - 5"
    end if

    if (.not.assert(qy%x, py%x)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_3 - 6"
    end if

    if (.not.assert(qy%y, py%y)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_3 - 7"
    end if

    if (.not.assert(qy%z, py%z)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_3 - 8"
    end if

    ! -----
    if (.not.assert(qz%w, pz%w)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_3 - 9"
    end if

    if (.not.assert(qz%x, pz%x)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_3 - 10"
    end if

    if (.not.assert(qz%y, pz%y)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_3 - 11"
    end if

    if (.not.assert(qz%z, pz%z)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_init_3 - 12"
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_abs() result(rst)
    logical :: rst
    type(quaternion) :: q
    real(real64) :: x(4), ans, qnorm

    ! Initialization
    rst = .true.
    call random_number(x)
    ans = norm2(x)
    q = quaternion(x)
    qnorm = abs(q)

    ! Test
    if (.not.assert(qnorm, ans)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_abs"
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_conjg() result(rst)
    logical :: rst
    type(quaternion) :: q, qc
    real(real64) :: x(4)

    ! Initialization
    rst = .true.
    call random_number(x)
    q = quaternion(x)
    qc = conjg(q)

    ! Test
    if (.not.assert(q%w, qc%w)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_conjg -1"
    end if
    if (.not.assert(q%x, -qc%x)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_conjg -2"
    end if
    if (.not.assert(q%y, -qc%y)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_conjg -3"
    end if
    if (.not.assert(q%z, -qc%z)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_conjg -4"
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_real() result(rst)
    logical :: rst
    type(quaternion) :: q
    real(real64) :: x(4)

    ! Initialization
    rst = .true.
    call random_number(x)
    q = quaternion(x)

    ! Test
    if (.not.assert(x(1), real(q))) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_real"
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_aimag() result(rst)
    logical :: rst
    type(quaternion) :: q
    real(real64) :: x(4), qimag(3)

    ! Initialization
    rst = .true.
    call random_number(x)
    q = quaternion(x)
    qimag = aimag(q)

    ! Test
    if (.not.assert(x(2:4), qimag)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_aimag"
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_multiply() result(rst)
    logical :: rst
    type(quaternion) :: q, qx, qy, q2
    real(real64) :: ax, ay, axisX(3), axisY(3), Rx(3,3), Ry(3,3), R(3,3), &
        rb(3), rq(3), rg(3), f

    ! Initialization
    rst = .true.
    call random_number(ax)
    call random_number(ay)
    call random_number(rb)
    call random_number(f)
    axisX = [1.0d0, 0.0d0, 0.0d0]
    axisY = [0.0d0, 1.0d0, 0.0d0]
    Rx = rotate_x(ax)
    Ry = rotate_y(ay)
    R = matmul(Rx, Ry)
    qx = quaternion(ax, axisX)
    qy = quaternion(ay, axisY)
    q = qx * qy
    call q%normalize()
    rg = matmul(R, rb)
    rq = aimag(q * rb * conjg(q))

    ! Test
    if (.not.assert(rg, rq)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_multiply -1"
    end if

    ! Scaling Test
    q2 = f * q
    if (.not.assert(q2%w, f * q%w)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_multiply -2"
    end if
    if (.not.assert(q2%x, f * q%x)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_multiply -3"
    end if
    if (.not.assert(q2%y, f * q%y)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_multiply -4"
    end if
    if (.not.assert(q2%z, f * q%z)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_multiply -5"
    end if

    ! Scaling Test 2
    q2 = q * f
    if (.not.assert(q2%w, f * q%w)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_multiply -6"
    end if
    if (.not.assert(q2%x, f * q%x)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_multiply -7"
    end if
    if (.not.assert(q2%y, f * q%y)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_multiply -8"
    end if
    if (.not.assert(q2%z, f * q%z)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_multiply -9"
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_add() result(rst)
    logical :: rst
    type(quaternion) :: q1, q2, q
    real(real64) :: x1(4), x2(4), x(4)

    ! Initialization
    rst = .true.
    call random_number(x1)
    call random_number(x2)
    x = x1 + x2
    q1 = quaternion(x1)
    q2 = quaternion(x2)
    q = q1 + q2

    ! Test
    if (.not.assert(x(1), q%w)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_add -1"
    end if

    if (.not.assert(x(2), q%x)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_add -2"
    end if

    if (.not.assert(x(3), q%y)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_add -3"
    end if

    if (.not.assert(x(4), q%z)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_add -4"
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_subtract() result(rst)
    logical :: rst
    type(quaternion) :: q1, q2, q
    real(real64) :: x1(4), x2(4), x(4)

    ! Initialization
    rst = .true.
    call random_number(x1)
    call random_number(x2)
    x = x1 - x2
    q1 = quaternion(x1)
    q2 = quaternion(x2)
    q = q1 - q2

    ! Test
    if (.not.assert(x(1), q%w)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_subtract -1"
    end if

    if (.not.assert(x(2), q%x)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_subtract -2"
    end if

    if (.not.assert(x(3), q%y)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_subtract -3"
    end if

    if (.not.assert(x(4), q%z)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_subtract -4"
    end if

    ! Negation Test
    x = -x
    q = -q
    if (.not.assert(x(1), q%w)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_subtract -5"
    end if

    if (.not.assert(x(2), q%x)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_subtract -6"
    end if

    if (.not.assert(x(3), q%y)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_subtract -7"
    end if

    if (.not.assert(x(4), q%z)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_subtract -8"
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_division() result(rst)
    logical :: rst
    type(quaternion) :: q1, q2, q, qa
    real(real64) :: x1(4), x2(4), f

    ! Initialization
    rst = .true.
    call random_number(x1)
    call random_number(x2)
    call random_number(f)
    q1 = quaternion(x1)
    q2 = quaternion(x2)
    q = q1 / q2
    qa = q1 * (conjg(q2) / (abs(q2)**2))

    ! Test
    if (.not.assert(q%w, qa%w)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_division -1"
    end if

    if (.not.assert(q%x, qa%x)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_division -2"
    end if

    if (.not.assert(q%y, qa%y)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_division -3"
    end if

    if (.not.assert(q%z, qa%z)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_division -4"
    end if

    ! Scalar Test
    qa = q / f
    if (.not.assert(qa%w, q%w / f)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_division -5"
    end if

    if (.not.assert(qa%x, q%x / f)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_division -6"
    end if

    if (.not.assert(qa%y, q%y / f)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_division -7"
    end if

    if (.not.assert(qa%z, q%z / f)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_division -8"
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_to_matrix() result(rst)
    logical :: rst
    type(quaternion) :: qx, qy, q
    real(real64) :: ax, ay, axisX(3), axisY(3), Rx(3,3), Ry(3,3), R(3,3), &
        Rq(3,3)
    
    ! Initialization
    rst = .true.
    call random_number(ax)
    call random_number(ay)
    axisX = [1.0d0, 0.0d0, 0.0d0]
    axisY = [0.0d0, 1.0d0, 0.0d0]
    Rx = rotate_x(ax)
    Ry = rotate_y(ay)
    R = matmul(Rx, Ry)
    qx = quaternion(ax, axisX)
    qy = quaternion(ay, axisY)
    q = qx * qy
    call q%normalize()
    Rq = q%to_matrix()

    ! Test
    if (.not.assert(R, Rq)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_to_matrix"
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_normalize() result(rst)
    logical :: rst
    type(quaternion) :: q, qn
    real(real64) :: x(4), xn(4)

    ! Initialization
    rst = .true.
    call random_number(x)
    xn = x / norm2(x)
    q = quaternion(x)
    qn = q
    call qn%normalize()

    ! Test
    if (.not.assert(xn(1), qn%w)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_normalize -1"
    end if

    if (.not.assert(xn(2), qn%x)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_normalize -2"
    end if

    if (.not.assert(xn(3), qn%y)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_normalize -3"
    end if

    if (.not.assert(xn(4), qn%z)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_normalize -4"
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_to_array() result(rst)
    logical :: rst
    type(quaternion) :: q
    real(real64) :: x(4)

    ! Initialization
    rst = .true.
    call random_number(x)
    q = quaternion(x)
    
    if (.not.assert(x, q%to_array())) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_to_array"
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_inverse() result(rst)
    logical :: rst
    type(quaternion) :: q, qi, qa
    real(real64) :: x(4)

    ! Initialization
    rst = .true.
    call random_number(x)
    q = quaternion(x)
    qi = inverse(q)
    qa = conjg(q) / (abs(q)**2)

    ! Test
    if (.not.assert(qi%w, qa%w)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_inverse -1"
    end if

    if (.not.assert(qi%x, qa%x)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_inverse -2"
    end if

    if (.not.assert(qi%y, qa%y)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_inverse -3"
    end if

    if (.not.assert(qi%z, qa%z)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_quaternion_inverse -4"
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_to_angle_axis() result(rst)
    logical :: rst
    real(real64) :: Rx(3,3), Ry(3,3), Rz(3,3), R(3,3), Ra(3,3), Rq(3,3), &
        axisX(3), axisY(3), axisZ(3), ax, ay, az, ar, aq, axisR(3), axisQ(3)
    type(quaternion) :: qx, qy, qz, q

    ! Initialization
    rst = .true.
    axisX = [1.0d0, 0.0d0, 0.0d0]
    axisY = [0.0d0, 1.0d0, 0.0d0]
    axisZ = [0.0d0, 0.0d0, 1.0d0]
    call random_number(ax)
    call random_number(ay)
    call random_number(az)
    qx = quaternion(ax, axisX)
    qy = quaternion(ay, axisY)
    qz = quaternion(az, axisZ)
    Rx = rotate_x(ax);
    Ry = rotate_y(ay);
    Rz = rotate_z(az);
    q = qx * qy * qz
    call q%normalize()
    R = matmul(Rx, matmul(Ry, Rz))

    ! Get the axis and angle information
    call to_angle_axis(R, ar, axisR)
    call q%to_angle_axis(aq, axisQ)

    ! Test
    if (.not.assert(ar, aq)) then
        print "(A)", "TEST FAILED: test_quaternion_to_angle_axis -1"
        rst = .false.
    end if
    if (.not.assert(axisR, axisQ)) then
        print "(A)", "TEST FAILED: test_quaternion_to_angle_axis -2"
        rst = .false.
    end if
end function

! ------------------------------------------------------------------------------
end module
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
module dynamics_rotation_tests
    use iso_fortran_env
    use fortran_test_helper
    use dynamics
    implicit none

contains
! ------------------------------------------------------------------------------
function test_axis_rotation_matrices() result(rst)
    !! Verifies the three Cartesian rotation matrices using known vector
    !! rotations and orthogonality.
    logical :: rst
    real(real64), parameter :: half_pi = 0.5d0 * acos(-1.0d0)
    real(real64), dimension(3,3) :: rx, ry, rz, identity

    rst = .true.
    identity = 0.0d0
    identity(1,1) = 1.0d0
    identity(2,2) = 1.0d0
    identity(3,3) = 1.0d0
    rx = rotate_x(half_pi)
    ry = rotate_y(half_pi)
    rz = rotate_z(half_pi)

    if (.not.assert(matmul(rx, [0.0d0, 1.0d0, 0.0d0]), &
        [0.0d0, 0.0d0, 1.0d0])) rst = .false.
    if (.not.assert(matmul(ry, [0.0d0, 0.0d0, 1.0d0]), &
        [1.0d0, 0.0d0, 0.0d0])) rst = .false.
    if (.not.assert(matmul(rz, [1.0d0, 0.0d0, 0.0d0]), &
        [0.0d0, 1.0d0, 0.0d0])) rst = .false.
    if (.not.assert(matmul(transpose(rx), rx), identity) .or. &
        .not.assert(matmul(transpose(ry), ry), identity) .or. &
        .not.assert(matmul(transpose(rz), rz), identity)) rst = .false.
    if (.not.rst) print "(A)", "TEST FAILED: test_axis_rotation_matrices"
end function

! ------------------------------------------------------------------------------
function test_homogeneous_rotation_matrices() result(rst)
    !! Verifies that each homogeneous rotation embeds its corresponding
    !! three-dimensional rotation without translation.
    logical :: rst
    real(real64), parameter :: angle = 0.37d0
    real(real64), dimension(4,4) :: hx, hy, hz
    real(real64), dimension(4) :: last_row

    rst = .true.
    last_row = [0.0d0, 0.0d0, 0.0d0, 1.0d0]
    hx = homogeneous_rotation_x(angle)
    hy = homogeneous_rotation_y(angle)
    hz = homogeneous_rotation_z(angle)
    if (.not.assert(hx(1:3,1:3), rotate_x(angle)) .or. &
        .not.assert(hy(1:3,1:3), rotate_y(angle)) .or. &
        .not.assert(hz(1:3,1:3), rotate_z(angle))) rst = .false.
    if (.not.assert(hx(1:3,4), [0.0d0, 0.0d0, 0.0d0]) .or. &
        .not.assert(hy(1:3,4), [0.0d0, 0.0d0, 0.0d0]) .or. &
        .not.assert(hz(1:3,4), [0.0d0, 0.0d0, 0.0d0])) rst = .false.
    if (.not.assert(hx(4,:), last_row) .or. &
        .not.assert(hy(4,:), last_row) .or. &
        .not.assert(hz(4,:), last_row)) rst = .false.
    if (.not.rst) print "(A)", &
        "TEST FAILED: test_homogeneous_rotation_matrices"
end function

! ------------------------------------------------------------------------------
function test_general_rotation_matrices() result(rst)
    !! Verifies both general frame-rotation overloads.
    logical :: rst
    real(real64), dimension(3,3) :: parent, child, actual

    rst = .true.
    child = matmul(rotate_z(0.4d0), rotate_y(-0.2d0))
    actual = rotate(child(:,1), child(:,2), child(:,3))
    if (.not.assert(actual, child)) rst = .false.

    parent = rotate_x(0.3d0)
    actual = rotate(child(:,1), child(:,2), child(:,3), &
        parent(:,1), parent(:,2), parent(:,3))
    if (.not.assert(actual, matmul(transpose(parent), child))) rst = .false.
    if (.not.rst) print "(A)", &
        "TEST FAILED: test_general_rotation_matrices"
end function

! ------------------------------------------------------------------------------
function test_translation_matrices() result(rst)
    !! Verifies both homogeneous translation overloads.
    logical :: rst
    real(real64), dimension(3) :: displacement
    real(real64), dimension(4,4) :: scalar_form, vector_form, identity

    rst = .true.
    displacement = [1.25d0, -2.5d0, 0.75d0]
    identity = 0.0d0
    identity(1,1) = 1.0d0
    identity(2,2) = 1.0d0
    identity(3,3) = 1.0d0
    identity(4,4) = 1.0d0
    scalar_form = translate(displacement(1), displacement(2), displacement(3))
    vector_form = translate(displacement)
    if (.not.assert(scalar_form, vector_form)) rst = .false.
    identity(1:3,4) = displacement
    if (.not.assert(vector_form, identity)) rst = .false.
    if (.not.rst) print "(A)", "TEST FAILED: test_translation_matrices"
end function

! ------------------------------------------------------------------------------
function test_rotation_angle_axis() result(rst)
    !! Verifies direct angle-axis extraction for a known non-singular rotation.
    logical :: rst
    real(real64), parameter :: expected_angle = 0.63d0
    real(real64) :: angle
    real(real64), dimension(3) :: axis

    rst = .true.
    call to_angle_axis(rotate_z(expected_angle), angle, axis)
    if (.not.assert(angle, expected_angle) .or. &
        .not.assert(axis, [0.0d0, 0.0d0, 1.0d0])) rst = .false.
    if (.not.rst) print "(A)", "TEST FAILED: test_rotation_angle_axis"
end function

! ------------------------------------------------------------------------------
function test_kinematic_transform_matrices() result(rst)
    !! Verifies velocity and acceleration transform matrices against their
    !! independent vector cross-product forms at an arbitrary point.
    logical :: rst
    real(real64), dimension(3) :: omega, alpha, velocity, acceleration, &
        origin, point, offset, expected_velocity, expected_acceleration
    real(real64), dimension(4) :: homogeneous_point, transformed
    real(real64), dimension(4,4) :: velocity_matrix, acceleration_matrix

    rst = .true.
    omega = [0.4d0, -0.2d0, 0.3d0]
    alpha = [-0.1d0, 0.5d0, 0.2d0]
    velocity = [1.0d0, -2.0d0, 0.5d0]
    acceleration = [-0.7d0, 0.8d0, 1.2d0]
    origin = [0.3d0, -0.6d0, 0.9d0]
    point = [-0.4d0, 1.1d0, 0.2d0]
    offset = point - origin
    homogeneous_point = [point, 1.0d0]

    velocity_matrix = velocity_transform(omega, velocity, origin)
    transformed = matmul(velocity_matrix, homogeneous_point)
    expected_velocity = velocity + cross_product(omega, offset)
    if (.not.assert(transformed(1:3), expected_velocity)) rst = .false.
    if (.not.assert(transformed(4), 0.0d0)) rst = .false.

    acceleration_matrix = acceleration_transform(alpha, omega, &
        acceleration, origin)
    transformed = matmul(acceleration_matrix, homogeneous_point)
    expected_acceleration = acceleration + cross_product(alpha, offset) + &
        cross_product(omega, cross_product(omega, offset))
    if (.not.assert(transformed(1:3), expected_acceleration)) rst = .false.
    if (.not.assert(transformed(4), 0.0d0)) rst = .false.
    if (.not.rst) print "(A)", &
        "TEST FAILED: test_kinematic_transform_matrices"
end function

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
function test_quaternion_exp() result(rst)
    logical :: rst
    type(quaternion) :: q, qexp
    real(real64) :: x(4), vmag, ans(4)

    ! Initialization
    rst = .true.
    call random_number(x)
    vmag = norm2(x(2:))
    ans(1) = exp(x(1)) * cos(vmag)
    ans(2:) = exp(x(1)) * x(2:) * sin(vmag) / vmag
    q = quaternion(x)
    qexp = exp(q)

    ! Test
    if (.not.assert(qexp%w, ans(1))) then
        print "(A)", "TEST FAILED: test_quaternion_exp -1"
        rst = .false.
    end if

    if (.not.assert(qexp%x, ans(2))) then
        print "(A)", "TEST FAILED: test_quaternion_exp -2"
        rst = .false.
    end if

    if (.not.assert(qexp%y, ans(3))) then
        print "(A)", "TEST FAILED: test_quaternion_exp -3"
        rst = .false.
    end if

    if (.not.assert(qexp%z, ans(4))) then
        print "(A)", "TEST FAILED: test_quaternion_exp -4"
        rst = .false.
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_log() result(rst)
    logical :: rst
    type(quaternion) :: q, qlog
    real(real64) :: x(4), qmag, vmag, ans(4)

    ! Initialization
    rst = .true.
    call random_number(x)
    vmag = norm2(x(2:))
    qmag = norm2(x)
    ans(1) = log(qmag)
    ans(2:) = x(2:) * acos(x(1) / qmag) / vmag
    q = quaternion(x)
    qlog = log(q)

    ! Test
    if (.not.assert(qlog%w, ans(1))) then
        print "(A)", "TEST FAILED: test_quaternion_log -1"
        rst = .false.
    end if

    if (.not.assert(qlog%x, ans(2))) then
        print "(A)", "TEST FAILED: test_quaternion_log -2"
        rst = .false.
    end if

    if (.not.assert(qlog%y, ans(3))) then
        print "(A)", "TEST FAILED: test_quaternion_log -3"
        rst = .false.
    end if

    if (.not.assert(qlog%z, ans(4))) then
        print "(A)", "TEST FAILED: test_quaternion_log -4"
        rst = .false.
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_pwr() result(rst)
    logical :: rst
    type(quaternion) :: q, qpow, q2
    real(real64) :: x(4), qmag, vmag, ans(4), phi
    real(real64), parameter :: e1 = 2.5d0
    integer(int32), parameter :: e2 = 2

    ! Initialization
    rst = .true.
    call random_number(x)
    vmag = norm2(x(2:))
    qmag = norm2(x)
    q = quaternion(x)
    phi = acos(q%w / qmag)

    ! Test 1
    ans(1) = (qmag**e1) * cos(e1 * phi)
    ans(2:) = (qmag**e1) * x(2:) * sin(e1 * phi) / vmag
    qpow = q**e1
    q2 = exp(log(q) * e1)
    if (.not.assert(qpow%w, ans(1))) then
        print "(A)", "TEST FAILED: test_quaternion_pwr -1"
        rst = .false.
    end if

    if (.not.assert(qpow%x, ans(2))) then
        print "(A)", "TEST FAILED: test_quaternion_pwr -2"
        rst = .false.
    end if

    if (.not.assert(qpow%y, ans(3))) then
        print "(A)", "TEST FAILED: test_quaternion_pwr -3"
        rst = .false.
    end if

    if (.not.assert(qpow%z, ans(4))) then
        print "(A)", "TEST FAILED: test_quaternion_pwr -4"
        rst = .false.
    end if

    ! Test 2
    ans(1) = (qmag**e2) * cos(e2 * phi)
    ans(2:) = (qmag**e2) * x(2:) * sin(e2 * phi) / vmag
    qpow = q**e2
    if (.not.assert(qpow%w, ans(1))) then
        print "(A)", "TEST FAILED: test_quaternion_pwr -5"
        rst = .false.
    end if

    if (.not.assert(qpow%x, ans(2))) then
        print "(A)", "TEST FAILED: test_quaternion_pwr -6"
        rst = .false.
    end if

    if (.not.assert(qpow%y, ans(3))) then
        print "(A)", "TEST FAILED: test_quaternion_pwr -7"
        rst = .false.
    end if

    if (.not.assert(qpow%z, ans(4))) then
        print "(A)", "TEST FAILED: test_quaternion_pwr -8"
        rst = .false.
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_dot_product() result(rst)
    logical :: rst
    type(quaternion) :: q1, q2
    real(real64) :: x1(4), x2(4), ans, dp

    ! Initialization
    rst = .true.
    call random_number(x1)
    call random_number(x2)
    ans = dot_product(x1, x2)
    q1 = quaternion(x1)
    q2 = quaternion(x2)
    dp = dot_product(q1, q2)

    ! Test
    if (.not.assert(dp, ans)) then
        print "(A)", "TEST FAILED: test_quaternion_dot_product"
        rst = .false.
    end if
end function

! ------------------------------------------------------------------------------
function test_quaternion_roll_pitch_yaw() result(rst)
    logical :: rst
    real(real64) :: Rx(3,3), Ry(3,3), Rz(3,3), R(3,3), ax, ay, az, rl, p, y
    type(quaternion) :: q

    ! Initialization
    rst = .true.
    call random_number(ax)
    call random_number(ay)
    call random_number(az)
    Rx = rotate_x(ax)
    Ry = rotate_y(ay)
    Rz = rotate_z(az)
    R = matmul(Rz, matmul(Ry, Rx))
    q = quaternion(R)
    call q%to_roll_pitch_yaw(rl, p, y)

    ! Test
    if (.not.assert(ax, rl)) then
        print "(A)", "TEST FAILED: test_quaternion_roll_pitch_yaw -1"
        rst = .false.
    end if

    if (.not.assert(ay, p)) then
        print "(A)", "TEST FAILED: test_quaternion_roll_pitch_yaw -2"
        rst = .false.
    end if

    if (.not.assert(az, y)) then
        print "(A)", "TEST FAILED: test_quaternion_roll_pitch_yaw -3"
        rst = .false.
    end if
end function

! ------------------------------------------------------------------------------
end module
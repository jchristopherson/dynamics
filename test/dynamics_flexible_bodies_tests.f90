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
module dynamics_flexible_bodies_tests
    use iso_fortran_env
    use ieee_arithmetic, only : ieee_is_finite
    use autodiff, only : hyperdual
    use dynamics
    implicit none

contains
! ------------------------------------------------------------------------------
function test_flexible_truss_2d() result(rst)
    !! Checks the 2D truss against closed-form axial force and stiffness.
    logical :: rst
    real(real64), parameter :: tol = 1.0d-9
    real(real64) :: l, e, fexp
    real(real64), dimension(4) :: u, f
    real(real64), dimension(4,4) :: k, kexp
    type(flexible_truss_2d) :: truss

    rst = .true.
    truss%reference_coordinates(:,1) = [0.0d0, 0.0d0]
    truss%reference_coordinates(:,2) = [1.0d0, 0.0d0]
    truss%area = 2.0d0
    truss%material = linear_elastic_1d(modulus = 1000.0d0)

    ! Stretched state: L = 1.01, E = (L^2 - L0^2) / (2 L0^2)
    u = [0.0d0, 0.0d0, 0.01d0, 0.0d0]
    call truss%internal_force(u, f, k)
    l = 1.01d0
    e = (l * l - 1.0d0) / 2.0d0
    fexp = 1000.0d0 * 2.0d0 * e * l
    if (abs(f(3) - fexp) > tol .or. abs(f(1) + fexp) > tol .or. &
        abs(f(2)) > tol .or. abs(f(4)) > tol) then
        print "(A)", "TEST FAILED: test_flexible_truss_2d - force"
        rst = .false.
        return
    end if
    if (maxval(abs(k - transpose(k))) > tol) then
        print "(A)", "TEST FAILED: test_flexible_truss_2d - symmetry"
        rst = .false.
        return
    end if

    ! Unstrained state: K = EA/L0 with the axial pattern [1 0 -1 0]
    u = 0.0d0
    call truss%internal_force(u, f, k)
    kexp = 0.0d0
    kexp(1,1) = 2000.0d0
    kexp(3,3) = 2000.0d0
    kexp(1,3) = -2000.0d0
    kexp(3,1) = -2000.0d0
    if (maxval(abs(f)) > tol .or. maxval(abs(k - kexp)) > tol) then
        print "(A)", "TEST FAILED: test_flexible_truss_2d - unstrained"
        rst = .false.
    end if
end function

! ------------------------------------------------------------------------------
function test_flexible_rigid_translation() result(rst)
    !! A rigid translation of a 3D truss must produce no internal force.
    logical :: rst
    real(real64), parameter :: tol = 1.0d-10
    real(real64), dimension(6) :: u, f
    real(real64), dimension(6,6) :: k
    type(flexible_truss_3d) :: truss

    rst = .true.
    truss%reference_coordinates(:,1) = [0.0d0, 0.0d0, 0.0d0]
    truss%reference_coordinates(:,2) = [1.0d0, 2.0d0, 2.0d0]
    truss%area = 1.5d0
    truss%material = linear_elastic_1d(modulus = 200.0d0)

    u = [1.0d0, 2.0d0, 3.0d0, 1.0d0, 2.0d0, 3.0d0]
    call truss%internal_force(u, f, k)
    if (maxval(abs(f)) > tol) then
        print "(A)", "TEST FAILED: test_flexible_rigid_translation"
        rst = .false.
    end if
end function

! ------------------------------------------------------------------------------
function test_flexible_truss_3d_fd() result(rst)
    !! Compares the AD tangent of a nonlinear 3D truss with central differences
    !! of its internal force.
    logical :: rst
    real(real64), parameter :: tol = 1.0d-5
    real(real64), dimension(6) :: u, f
    real(real64), dimension(6,6) :: k, kfd
    type(flexible_truss_3d) :: truss

    rst = .true.
    truss%reference_coordinates(:,1) = [0.0d0, 0.0d0, 0.0d0]
    truss%reference_coordinates(:,2) = [1.0d0, 2.0d0, 2.0d0]
    truss%area = 1.5d0
    truss%material = polynomial_elastic_1d(modulus = 200.0d0, &
        quadratic = 5.0d0, quartic = 2.0d0)

    u = [0.1d0, -0.05d0, 0.2d0, 0.3d0, 0.2d0, -0.1d0]
    call truss%internal_force(u, f, k)
    call fd_stiffness(truss, u, kfd)
    if (maxval(abs(k - kfd)) > tol * max(1.0d0, maxval(abs(k)))) then
        print "(A)", "TEST FAILED: test_flexible_truss_3d_fd"
        rst = .false.
    end if
end function

! ------------------------------------------------------------------------------
function test_flexible_spring() result(rst)
    !! Checks a nonlinear 2D spring against closed-form force and stiffness,
    !! then checks a 3D spring against finite differences.
    logical :: rst
    real(real64), parameter :: tol = 1.0d-9
    real(real64) :: delta, fexp, kexp
    real(real64), dimension(4) :: u, f
    real(real64), dimension(4,4) :: k
    real(real64), dimension(6) :: u3, f3
    real(real64), dimension(6,6) :: k3, kfd
    type(flexible_spring_2d) :: spring
    type(flexible_spring_3d) :: spring3

    rst = .true.
    spring%reference_coordinates(:,1) = [0.0d0, 0.0d0]
    spring%reference_coordinates(:,2) = [1.0d0, 0.0d0]
    spring%stiffness = 10.0d0
    spring%cubic_stiffness = 5.0d0

    u = [0.0d0, 0.0d0, 0.2d0, 0.0d0]
    call spring%internal_force(u, f, k)
    delta = 0.2d0
    fexp = 10.0d0 * delta + 5.0d0 * delta**3
    kexp = 10.0d0 + 3.0d0 * 5.0d0 * delta * delta
    if (abs(f(3) - fexp) > tol .or. abs(f(1) + fexp) > tol .or. &
        abs(k(3,3) - kexp) > tol .or. abs(k(1,3) + kexp) > tol) then
        print "(A)", "TEST FAILED: test_flexible_spring - 2D"
        rst = .false.
        return
    end if

    spring3%reference_coordinates(:,1) = [0.0d0, 0.0d0, 0.0d0]
    spring3%reference_coordinates(:,2) = [1.0d0, 2.0d0, 2.0d0]
    spring3%stiffness = 3.0d0
    spring3%cubic_stiffness = 1.5d0
    u3 = [0.05d0, -0.1d0, 0.02d0, 0.2d0, 0.1d0, -0.3d0]
    call spring3%internal_force(u3, f3, k3)
    call fd_stiffness(spring3, u3, kfd)
    if (maxval(abs(k3 - kfd)) > 1.0d-5 * max(1.0d0, maxval(abs(k3)))) then
        print "(A)", "TEST FAILED: test_flexible_spring - 3D"
        rst = .false.
    end if
end function

! ------------------------------------------------------------------------------
function test_flexible_damper() result(rst)
    !! Checks the damping matrix of 2D and 3D viscous dampers at the reference
    !! configuration, where C = c b b^T with b = [-n, n] and n the unit axis.
    logical :: rst
    real(real64), parameter :: tol = 1.0d-10
    real(real64), dimension(4) :: u2
    real(real64), dimension(6) :: u3
    real(real64), dimension(4) :: b2
    real(real64), dimension(6) :: b3
    real(real64), dimension(4,4) :: c2, c2exp
    real(real64), dimension(6,6) :: c3, c3exp
    integer(int32) :: row, col
    type(flexible_damper_2d) :: damper2
    type(flexible_damper_3d) :: damper3

    rst = .true.
    damper2%reference_coordinates(:,1) = [0.0d0, 0.0d0]
    damper2%reference_coordinates(:,2) = [3.0d0, 4.0d0]
    damper2%damping = 2.0d0
    u2 = 0.0d0
    call damper2%damping_matrix(u2, c2)
    b2 = [-0.6d0, -0.8d0, 0.6d0, 0.8d0]
    do col = 1, 4
        do row = 1, 4
            c2exp(row,col) = 2.0d0 * b2(row) * b2(col)
        end do
    end do
    if (maxval(abs(c2 - c2exp)) > tol) then
        print "(A)", "TEST FAILED: test_flexible_damper - 2D"
        rst = .false.
    end if

    damper3%reference_coordinates(:,1) = [0.0d0, 0.0d0, 0.0d0]
    damper3%reference_coordinates(:,2) = [1.0d0, 2.0d0, 2.0d0]
    damper3%damping = 0.5d0
    u3 = 0.0d0
    call damper3%damping_matrix(u3, c3)
    b3 = [-1.0d0/3.0d0, -2.0d0/3.0d0, -2.0d0/3.0d0, &
        1.0d0/3.0d0, 2.0d0/3.0d0, 2.0d0/3.0d0]
    do col = 1, 6
        do row = 1, 6
            c3exp(row,col) = 0.5d0 * b3(row) * b3(col)
        end do
    end do
    if (maxval(abs(c3 - c3exp)) > tol) then
        print "(A)", "TEST FAILED: test_flexible_damper - 3D"
        rst = .false.
    end if
end function

! ------------------------------------------------------------------------------
function test_flexible_beam_linear_limits() result(rst)
    !! Compares unloaded tangents against the classical planar bending and
    !! axial stiffnesses and checks spatial torsional stiffness independently.
    logical :: rst
    type(flexible_euler_bernoulli_2d) :: eb
    type(flexible_timoshenko_2d) :: tm
    type(flexible_euler_bernoulli_3d) :: eb3
    type(flexible_timoshenko_3d) :: tm3
    real(real64) :: f(6), k(6,6), f3(12), k3(12,12), expected(2,2)
    integer, parameter :: tip(2) = [5, 6]
    real(real64), parameter :: l = 2.0d0, ei = 3.0d0, ga = 7.0d0, gj = 11.0d0
    real(real64), parameter :: tol = 1.0d-10, twist = 0.6d0
    real(real64) :: u3(12)

    eb%reference_coordinates(:,2) = [l, 0.0d0]
    eb%bending_stiffness = ei
    eb%area = 2.0d0
    eb%material = linear_elastic_1d(5.0d0)
    tm%reference_coordinates = eb%reference_coordinates
    tm%bending_stiffness = ei
    tm%shear_stiffness = ga
    tm%area = eb%area
    tm%material = eb%material

    call eb%internal_force([0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0], f, k)
    expected = ei * reshape([12/l**3, -6/l**2, -6/l**2, 4/l], [2,2])
    rst = all(ieee_is_finite(k)) .and. maxval(abs(k(tip,tip) - expected)) < tol .and. &
        abs(k(4,4) - 10/l) < tol .and. maxval(abs(f)) < tol
    if (.not. rst) then
        print *, "TEST FAILED: Euler-Bernoulli planar linear stiffness"
        return
    end if

    call tm%internal_force([0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0], f, k)
    expected = reshape([ga/l, -ga/2, -ga/2, ei/l + ga*l/4], [2,2])
    rst = all(ieee_is_finite(k)) .and. maxval(abs(k(tip,tip) - expected)) < tol .and. &
        abs(k(4,4) - 10/l) < tol .and. maxval(abs(f)) < tol
    if (.not. rst) then
        print *, "TEST FAILED: Timoshenko planar linear stiffness"
        return
    end if

    eb3%reference_coordinates(:,2) = [l, 0.0d0, 0.0d0]
    eb3%torsional_stiffness = gj
    eb3%material = eb%material
    tm3%reference_coordinates = eb3%reference_coordinates
    tm3%torsional_stiffness = gj
    tm3%material = eb%material
    u3 = 0.0d0
    call eb3%internal_force(u3, f3, k3)
    rst = all(ieee_is_finite(k3)) .and. abs(k3(10,10) - gj/l) < tol
    u3(10) = twist
    call tm3%internal_force(u3, f3, k3)
    rst = rst .and. all(ieee_is_finite(k3)) .and. abs(f3(10) - gj*twist/l) < tol .and. &
        abs(k3(10,10) - gj/l) < tol
    if (.not. rst) print *, "TEST FAILED: spatial beam torsion"
end function

function test_flexible_beam_objectivity() result(rst)
    !! Verifies that a finite rigid motion gives no elastic force and that
    !! superposing that motion on a strained state preserves its energy.
    logical :: rst
    type(flexible_euler_bernoulli_3d) :: eb
    type(flexible_timoshenko_3d) :: tm
    real(real64), parameter :: axis(3) = [2.0d0, -1.0d0, 2.0d0]/3.0d0
    real(real64), parameter :: shift(3) = [0.2d0, -0.5d0, 0.4d0], theta = 1.1d0
    real(real64) :: u(12), rigid(12), moved(12), p(3), r(3), unit_axis(3), angle
    type(quaternion) :: superposed
    integer :: i, j

    eb%reference_coordinates(:,1) = [0.2d0, 0.1d0, -0.4d0]
    eb%reference_coordinates(:,2) = [1.5d0, 1.0d0, 0.6d0]
    eb%bending_stiffness = [2.0d0, 5.0d0]
    eb%material = polynomial_elastic_1d(20.0d0, 9.0d0, 3.0d0)
    tm%reference_coordinates = eb%reference_coordinates
    tm%bending_stiffness = eb%bending_stiffness
    tm%shear_stiffness = [3.0d0, 7.0d0]
    tm%material = eb%material
    u = [0.01d0, 0.04d0, -0.03d0, -0.15d0, 0.1d0, 0.05d0, &
        0.2d0, -0.08d0, 0.1d0, 0.2d0, 0.15d0, -0.25d0]
    do i = 1, 2
        j = 6*(i-1)
        p = eb%reference_coordinates(:,i)
        rigid(j+1:j+3) = rotate_test_point(p, axis, theta) + shift - p
        rigid(j+4:j+6) = theta*axis
        moved(j+1:j+3) = rotate_test_point(p + u(j+1:j+3), axis, theta) + shift - p
        r = u(j+4:j+6)
        superposed = quaternion(theta, axis) * quaternion(norm2(r), r/norm2(r))
        call superposed%to_angle_axis(angle, unit_axis)
        moved(j+4:j+6) = angle*unit_axis
    end do
    rst = check_beam_objectivity(eb, u, rigid, moved) .and. check_beam_objectivity(tm, u, rigid, moved)
    if (.not. rst) print *, "TEST FAILED: flexible beam finite-motion objectivity"
end function

function rotate_test_point(p, axis, angle) result(q)
    !! Uses Rodrigues' formula to construct the independent rigid-motion test.
    real(real64), intent(in) :: p(3), axis(3), angle
    real(real64) :: q(3)

    q = cos(angle)*p + sin(angle)*cross_product(axis, p) + (1-cos(angle))*dot_product(axis,p)*axis
end function

function check_beam_objectivity(el, u, rigid, moved) result(rst)
    !! Checks energy invariance and compares the tangent at a strained state
    !! with central differences of the internal force.
    class(finite_element), intent(in) :: el
    real(real64), intent(in) :: u(12), rigid(12), moved(12)
    logical :: rst
    real(real64) :: f(12), k(12,12), kfd(12,12)
    type(hyperdual) :: seed(12), w, moved_w, rigid_w
    integer :: i

    do i = 1, 12
        seed(i) = hyperdual(u(i), 0.0d0, 0.0d0, 0.0d0)
    end do
    w = el%strain_energy(seed)
    seed%x = moved
    moved_w = el%strain_energy(seed)
    seed%x = rigid
    rigid_w = el%strain_energy(seed)
    call el%internal_force(rigid, f, k)
    rst = el%get_dof_count() == 12 .and. abs(rigid_w%x) < 1.0d-12 .and. &
        maxval(abs(f)) < 1.0d-10 .and. w%x > 1.0d-5 .and. abs(w%x - moved_w%x) < 1.0d-11
    call el%internal_force(moved, f, k)
    call fd_stiffness(el, moved, kfd)
    rst = rst .and. all(ieee_is_finite(k)) .and. all(ieee_is_finite(kfd)) .and. &
        maxval(abs(k - kfd)) < 1.0d-6 * max(1.0d0, maxval(abs(k)))
end function

! ------------------------------------------------------------------------------
subroutine fd_stiffness(el, u, kfd)
    !! Estimates the tangent stiffness by central differences of the internal
    !! force of a finite element.
    class(finite_element), intent(in) :: el
    real(real64), intent(in), dimension(:) :: u
    real(real64), intent(out), dimension(:,:) :: kfd

    real(real64), parameter :: h = 1.0d-6
    integer(int32) :: j, n
    real(real64), allocatable, dimension(:) :: up, fp, fm
    real(real64), allocatable, dimension(:,:) :: kdummy

    n = size(u)
    allocate(up(n), fp(n), fm(n), kdummy(n,n))
    do j = 1, n
        up = u
        up(j) = u(j) + h
        call el%internal_force(up, fp, kdummy)
        up = u
        up(j) = u(j) - h
        call el%internal_force(up, fm, kdummy)
        kfd(:,j) = (fp - fm) / (2.0d0 * h)
    end do
end subroutine

! ------------------------------------------------------------------------------
end module

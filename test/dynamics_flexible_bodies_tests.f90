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

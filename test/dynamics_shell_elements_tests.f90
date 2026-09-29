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
module dynamics_shell_elements_tests
    use iso_fortran_env
    use dynamics
    implicit none

    ! Patch-test displacement field parameters: membrane strains (a, b, c, d)
    ! and Kirchhoff curvature parameters (k1, k2, k3).
    real(real64), parameter, private :: pa = 1.0d-3
    real(real64), parameter, private :: pb = 2.0d-4
    real(real64), parameter, private :: pc = -3.0d-4
    real(real64), parameter, private :: pd = 5.0d-4
    real(real64), parameter, private :: pk1 = 2.0d-2
    real(real64), parameter, private :: pk2 = -1.0d-2
    real(real64), parameter, private :: pk3 = 1.5d-2

contains
! ******************************************************************************
! HELPERS
! ------------------------------------------------------------------------------
pure function rigid_body_vectors(p) result(rst)
    ! Builds the six rigid-body displacement vectors (three translations and
    ! three rotations about the global origin) for nodes located at p(3,n).
    real(real64), intent(in), dimension(:,:) :: p
    real(real64), allocatable, dimension(:,:) :: rst

    integer(int32) :: i, j, c
    real(real64) :: w(3)

    allocate(rst(6 * size(p, 2), 6), source = 0.0d0)
    do i = 1, size(p, 2)
        c = 6 * (i - 1)
        do j = 1, 3
            rst(c+j,j) = 1.0d0
            w = 0.0d0
            w(j) = 1.0d0
            rst(c+1:c+3,3+j) = cross_product(w, p(:,i))
            rst(c+4:c+6,3+j) = w
        end do
    end do
end function

! ------------------------------------------------------------------------------
pure function patch_displacement(x, y) result(rst)
    ! Nodal [u, v, w, theta_x, theta_y, theta_z] for a planar patch field
    ! with constant membrane strain and constant Kirchhoff curvature.
    real(real64), intent(in) :: x, y
    real(real64) :: rst(6)
    rst(1) = pa * x + pb * y
    rst(2) = pc * x + pd * y
    rst(3) = 0.5d0 * pk1 * x**2 + 0.5d0 * pk2 * y**2 + pk3 * x * y
    rst(4) = pk2 * y + pk3 * x
    rst(5) = -(pk1 * x + pk3 * y)
    rst(6) = 0.5d0 * (pc - pb)
end function

! ------------------------------------------------------------------------------
pure function patch_strain() result(rst)
    ! The exact generalized strain vector (first six components) produced by
    ! patch_displacement.
    real(real64) :: rst(6)
    rst = [pa, pd, pb + pc, -pk1, -pk2, -2.0d0 * pk3]
end function

! ------------------------------------------------------------------------------
subroutine build_strip(nx, ny, lx, ly, plane_xz, nodes)
    ! Builds a regular (nx+1)-by-(ny+1) grid of 6-DOF nodes lying in either
    ! the global x-y plane or the global x-z plane.
    integer(int32), intent(in) :: nx, ny
    real(real64), intent(in) :: lx, ly
    logical, intent(in) :: plane_xz
    type(node), allocatable, intent(out), dimension(:) :: nodes

    integer(int32) :: i, j, k
    real(real64) :: x, y

    allocate(nodes((nx + 1) * (ny + 1)))
    do j = 1, ny + 1
        do i = 1, nx + 1
            k = (j - 1) * (nx + 1) + i
            x = lx * (i - 1) / nx
            y = ly * (j - 1) / ny
            if (plane_xz) then
                nodes(k) = node(k, 6, x, 0.0d0, y)
            else
                nodes(k) = node(k, 6, x, y, 0.0d0)
            end if
        end do
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine clamp_and_load(nx, ny, dir, bc, f)
    ! Clamps all DOF along x = 0 and applies a unit total tip load in the
    ! global direction dir, distributed consistently along the tip edge.
    integer(int32), intent(in) :: nx, ny, dir
    integer(int32), allocatable, intent(out), dimension(:) :: bc
    real(real64), intent(inout), dimension(:) :: f

    integer(int32) :: j, m, kroot, ktip
    real(real64) :: wt

    allocate(bc(6 * (ny + 1)))
    f = 0.0d0
    do j = 1, ny + 1
        kroot = (j - 1) * (nx + 1) + 1
        ktip = kroot + nx
        do m = 1, 6
            bc(6 * (j - 1) + m) = 6 * (kroot - 1) + m
        end do
        wt = 1.0d0 / ny
        if (j == 1 .or. j == ny + 1) wt = 0.5d0 * wt
        f(6 * (ktip - 1) + dir) = wt
    end do
end subroutine

! ------------------------------------------------------------------------------
function cantilever_tip_deflection(elements, nodes, nx, ny, dir) result(rst)
    ! Solves the clamped strip under a unit tip load and returns the mean
    ! tip deflection in the global direction dir.
    class(element), intent(in), dimension(:) :: elements
    type(node), intent(in), dimension(:) :: nodes
    integer(int32), intent(in) :: nx, ny, dir
    real(real64) :: rst

    integer(int32) :: j, gdof, ktip
    integer(int32), allocatable, dimension(:) :: bc
    real(real64), allocatable, dimension(:) :: f, fr, ur, u
    real(real64), allocatable, dimension(:,:) :: k, kr

    gdof = 6 * size(nodes)
    call assemble_static_system(gdof, elements, nodes, k)
    allocate(f(gdof))
    call clamp_and_load(nx, ny, dir, bc, f)
    kr = apply_boundary_conditions(bc, k)
    fr = apply_boundary_conditions(bc, f)
    ur = solve_static_system(kr, fr)
    u = restore_constrained_values(bc, ur)

    rst = 0.0d0
    do j = 1, ny + 1
        ktip = (j - 1) * (nx + 1) + nx + 1
        rst = rst + u(6 * (ktip - 1) + dir)
    end do
    rst = rst / (ny + 1)
end function

! ------------------------------------------------------------------------------
function cantilever_first_frequency(elements, nodes, nx, ny) result(rst)
    ! Computes the first natural frequency (rad/s) of the clamped strip.
    class(element), intent(in), dimension(:) :: elements
    type(node), intent(in), dimension(:) :: nodes
    integer(int32), intent(in) :: nx, ny
    real(real64) :: rst

    integer(int32) :: gdof
    integer(int32), allocatable, dimension(:) :: bc
    real(real64), allocatable, dimension(:) :: f, freqs
    real(real64), allocatable, dimension(:,:) :: m, k, mr, kr

    gdof = 6 * size(nodes)
    call assemble_dynamic_system(gdof, elements, nodes, m, k)
    allocate(f(gdof))
    call clamp_and_load(nx, ny, 3, bc, f)
    mr = apply_boundary_conditions(bc, m)
    kr = apply_boundary_conditions(bc, k)
    call modal_response(mr, kr, freqs)
    rst = freqs(1)
end function

! ******************************************************************************
! TESTS
! ------------------------------------------------------------------------------
function test_shell_rigid_body_modes() result(rst)
    logical :: rst
    integer(int32) :: j
    real(real64), parameter :: tol = 1.0d-9
    real(real64) :: kscale, energy
    real(real64), allocatable, dimension(:) :: u
    real(real64), allocatable, dimension(:,:) :: k, rb, p
    type(material) :: mat
    type(node) :: n1, n2, n3, n4
    type(triangular_shell_element) :: tri
    type(rectangular_shell_element) :: quad

    rst = .true.
    mat = material(2.0d11, 0.3d0, 7.85d3)

    ! Arbitrarily oriented triangle
    n1 = node(1, 6, 0.1d0, 0.2d0, 0.3d0)
    n2 = node(2, 6, 1.3d0, 0.4d0, 0.9d0)
    n3 = node(3, 6, 0.4d0, 1.1d0, 1.2d0)
    tri = triangular_shell_element(mat, 0.02d0, n1, n2, n3)
    k = tri%stiffness_matrix()
    p = reshape([0.1d0, 0.2d0, 0.3d0, 1.3d0, 0.4d0, 0.9d0, &
        0.4d0, 1.1d0, 1.2d0], [3, 3])
    rb = rigid_body_vectors(p)
    kscale = maxval(abs(k))
    if (size(k, 1) /= 18 .or. .not.is_symmetric(k)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_shell_rigid_body_modes - tri symmetry"
    end if
    do j = 1, 6
        if (norm2(matmul(k, rb(:,j))) > tol * kscale * norm2(rb(:,j))) then
            rst = .false.
            print "(A, I0)", &
                "TEST FAILED: test_shell_rigid_body_modes - tri mode ", j
        end if
    end do

    ! A deformation mode must store strain energy
    allocate(u(18), source = 0.0d0)
    u(3) = 1.0d0
    energy = dot_product(u, matmul(k, u))
    if (energy <= 0.0d0) then
        rst = .false.
        print "(A)", "TEST FAILED: test_shell_rigid_body_modes - tri energy"
    end if

    ! Arbitrarily oriented, planar, non-rectangular quadrilateral lying on
    ! the plane z = 0.5 x + 0.25 y + 0.1
    p = reshape([0.0d0, 0.0d0, 0.1d0, 2.0d0, 0.2d0, 1.15d0, &
        1.8d0, 1.5d0, 1.375d0, 0.3d0, 1.2d0, 0.55d0], [3, 4])
    n1 = node(1, 6, p(1,1), p(2,1), p(3,1))
    n2 = node(2, 6, p(1,2), p(2,2), p(3,2))
    n3 = node(3, 6, p(1,3), p(2,3), p(3,3))
    n4 = node(4, 6, p(1,4), p(2,4), p(3,4))
    quad = rectangular_shell_element(mat, 0.05d0, n1, n2, n3, n4)
    k = quad%stiffness_matrix()
    rb = rigid_body_vectors(p)
    kscale = maxval(abs(k))
    if (size(k, 1) /= 24 .or. .not.is_symmetric(k)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_shell_rigid_body_modes - quad symmetry"
    end if
    do j = 1, 6
        if (norm2(matmul(k, rb(:,j))) > tol * kscale * norm2(rb(:,j))) then
            rst = .false.
            print "(A, I0)", &
                "TEST FAILED: test_shell_rigid_body_modes - quad mode ", j
        end if
    end do
    deallocate(u)
    allocate(u(24), source = 0.0d0)
    u(3) = 1.0d0
    energy = dot_product(u, matmul(k, u))
    if (energy <= 0.0d0) then
        rst = .false.
        print "(A)", "TEST FAILED: test_shell_rigid_body_modes - quad energy"
    end if
end function

! ------------------------------------------------------------------------------
function test_shell_patch_tests() result(rst)
    logical :: rst
    integer(int32) :: i, j
    real(real64), parameter :: tol = 1.0d-10
    real(real64), parameter :: e = 7.0d10
    real(real64), parameter :: nu = 0.3d0
    real(real64), parameter :: t = 0.01d0
    real(real64) :: xy(2,6), expected(6), sexp(6), c0
    real(real64), allocatable, dimension(:) :: u, strain, stress
    real(real64), allocatable, dimension(:,:) :: pts, avg
    type(material) :: mat
    type(node), dimension(6) :: nodes
    type(triangular_shell_element) :: tri
    type(rectangular_shell_element), dimension(2) :: quads

    rst = .true.
    mat = material(e, nu, 2.7d3)
    expected = patch_strain()

    ! Expected stress resultants
    c0 = e / (1.0d0 - nu**2)
    sexp(1) = t * c0 * (expected(1) + nu * expected(2))
    sexp(2) = t * c0 * (nu * expected(1) + expected(2))
    sexp(3) = t * c0 * 0.5d0 * (1.0d0 - nu) * expected(3)
    sexp(4) = (t**3 / 12.0d0) * c0 * (expected(4) + nu * expected(5))
    sexp(5) = (t**3 / 12.0d0) * c0 * (nu * expected(4) + expected(5))
    sexp(6) = (t**3 / 12.0d0) * c0 * 0.5d0 * (1.0d0 - nu) * expected(6)

    ! Distorted triangle: DKT reproduces constant curvature exactly and has
    ! no transverse shear.
    xy(:,1:3) = reshape([0.0d0, 0.0d0, 2.0d0, 0.0d0, 0.5d0, 1.5d0], [2, 3])
    do i = 1, 3
        nodes(i) = node(i, 6, xy(1,i), xy(2,i), 0.0d0)
    end do
    tri = triangular_shell_element(mat, t, nodes(1), nodes(2), nodes(3))
    allocate(u(18))
    do i = 1, 3
        u(6*i-5:6*i) = patch_displacement(xy(1,i), xy(2,i))
    end do
    pts = reshape([1.0d0 / 3.0d0, 1.0d0 / 3.0d0, 0.1d0, 0.7d0, &
        0.0d0, 0.0d0, 0.0d0, 1.0d0], [2, 4])
    do j = 1, size(pts, 2)
        strain = tri%strain(u, pts(:,j))
        stress = tri%stress(u, pts(:,j))
        if (size(strain) /= 8 .or. &
            maxval(abs(strain(1:6) - expected)) > tol .or. &
            maxval(abs(strain(7:8))) > tol .or. &
            maxval(abs(stress(1:6) - sexp) / max(abs(sexp), 1.0d0)) > &
            1.0d-8) then
            rst = .false.
            print "(A, I0)", "TEST FAILED: test_shell_patch_tests - tri point ", j
        end if
    end do

    ! Two distorted quadrilaterals sharing an edge
    xy = reshape([0.0d0, 0.0d0, 1.0d0, 0.0d0, 2.0d0, 0.0d0, &
        -0.1d0, 1.0d0, 1.2d0, 1.1d0, 2.1d0, 1.3d0], [2, 6])
    do i = 1, 6
        nodes(i) = node(i, 6, xy(1,i), xy(2,i), 0.0d0)
    end do
    quads(1) = rectangular_shell_element(mat, t, nodes(1), nodes(2), &
        nodes(5), nodes(4))
    quads(2) = rectangular_shell_element(mat, t, nodes(2), nodes(3), &
        nodes(6), nodes(5))
    deallocate(u)
    allocate(u(36))
    do i = 1, 6
        u(6*i-5:6*i) = patch_displacement(xy(1,i), xy(2,i))
    end do

    ! Pointwise checks on the first quadrilateral (membrane and curvature
    ! fields are linear in the rotations and are reproduced exactly)
    pts = reshape([0.0d0, 0.0d0, 0.5d0, -0.3d0, 1.0d0, 1.0d0], [2, 3])
    do j = 1, size(pts, 2)
        strain = quads(1)%strain([u(1:12), u(25:30), u(19:24)], pts(:,j))
        stress = quads(1)%stress([u(1:12), u(25:30), u(19:24)], pts(:,j))
        if (size(strain) /= 8 .or. &
            maxval(abs(strain(1:6) - expected)) > tol .or. &
            maxval(abs(stress(1:6) - sexp) / max(abs(sexp), 1.0d0)) > &
            1.0d-8) then
            rst = .false.
            print "(A, I0)", "TEST FAILED: test_shell_patch_tests - quad point ", j
        end if
    end do

    ! Nodal averaging over the patch
    avg = nodally_averaged_strain(quads, nodes, u)
    if (size(avg, 1) /= 8 .or. size(avg, 2) /= 6) then
        rst = .false.
        print "(A)", "TEST FAILED: test_shell_patch_tests - averaged size"
    else
        do i = 1, 6
            if (maxval(abs(avg(1:6,i) - expected)) > tol) then
                rst = .false.
                print "(A, I0)", &
                    "TEST FAILED: test_shell_patch_tests - averaged node ", i
            end if
        end do
    end if
end function

! ------------------------------------------------------------------------------
function test_shell_mass_and_load() result(rst)
    logical :: rst
    integer(int32) :: i
    real(real64), parameter :: tol = 1.0d-10
    real(real64), parameter :: rho = 7.85d3
    real(real64), parameter :: t = 0.02d0
    real(real64) :: area, axis(3), q(3), ftot(3)
    real(real64), allocatable, dimension(:) :: u, f
    real(real64), allocatable, dimension(:,:) :: m, p
    type(material) :: mat
    type(node) :: n1, n2, n3, n4
    type(triangular_shell_element) :: tri
    type(rectangular_shell_element) :: quad

    rst = .true.
    mat = material(2.0d11, 0.3d0, rho)
    axis = [1.0d0, 2.0d0, 2.0d0] / 3.0d0
    q = [1.0d0, -2.0d0, 3.0d0]

    ! Triangle
    p = reshape([0.1d0, 0.2d0, 0.3d0, 1.3d0, 0.4d0, 0.9d0, &
        0.4d0, 1.1d0, 1.2d0], [3, 3])
    n1 = node(1, 6, p(1,1), p(2,1), p(3,1))
    n2 = node(2, 6, p(1,2), p(2,2), p(3,2))
    n3 = node(3, 6, p(1,3), p(2,3), p(3,3))
    tri = triangular_shell_element(mat, t, n1, n2, n3)
    area = 0.5d0 * norm2(cross_product(p(:,2) - p(:,1), p(:,3) - p(:,1)))
    m = tri%mass_matrix()
    if (abs(tri%area() - area) > tol * area .or. .not.is_symmetric(m)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_shell_mass_and_load - tri area/symmetry"
    end if

    ! Rigid translation carries the total mass; rigid rotation of the
    ! rotational DOF carries the rotary inertia
    allocate(u(18), source = 0.0d0)
    do i = 1, 3
        u(6*i-5:6*i-3) = axis
    end do
    if (abs(dot_product(u, matmul(m, u)) - rho * t * area) > &
        tol * rho * t * area) then
        rst = .false.
        print "(A)", "TEST FAILED: test_shell_mass_and_load - tri mass"
    end if
    u = 0.0d0
    do i = 1, 3
        u(6*i-2:6*i) = axis
    end do
    if (abs(dot_product(u, matmul(m, u)) - rho * t**3 * area / 12.0d0) > &
        tol * rho * t**3 * area / 12.0d0) then
        rst = .false.
        print "(A)", "TEST FAILED: test_shell_mass_and_load - tri inertia"
    end if

    ! Uniform traction is shared equally among the triangle nodes
    f = tri%external_force_vector(q)
    do i = 1, 3
        if (maxval(abs(f(6*i-5:6*i-3) - area * q / 3.0d0)) > tol .or. &
            maxval(abs(f(6*i-2:6*i))) > 0.0d0) then
            rst = .false.
            print "(A, I0)", "TEST FAILED: test_shell_mass_and_load - tri load ", i
        end if
    end do

    ! Planar quadrilateral
    p = reshape([0.0d0, 0.0d0, 0.1d0, 2.0d0, 0.2d0, 1.15d0, &
        1.8d0, 1.5d0, 1.375d0, 0.3d0, 1.2d0, 0.55d0], [3, 4])
    n1 = node(1, 6, p(1,1), p(2,1), p(3,1))
    n2 = node(2, 6, p(1,2), p(2,2), p(3,2))
    n3 = node(3, 6, p(1,3), p(2,3), p(3,3))
    n4 = node(4, 6, p(1,4), p(2,4), p(3,4))
    quad = rectangular_shell_element(mat, t, n1, n2, n3, n4)
    area = 0.5d0 * norm2(cross_product(p(:,3) - p(:,1), p(:,4) - p(:,2)))
    m = quad%mass_matrix()
    if (abs(quad%area() - area) > tol * area .or. .not.is_symmetric(m)) then
        rst = .false.
        print "(A)", &
            "TEST FAILED: test_shell_mass_and_load - quad area/symmetry"
    end if
    deallocate(u)
    allocate(u(24), source = 0.0d0)
    do i = 1, 4
        u(6*i-5:6*i-3) = axis
    end do
    if (abs(dot_product(u, matmul(m, u)) - rho * t * area) > &
        tol * rho * t * area) then
        rst = .false.
        print "(A)", "TEST FAILED: test_shell_mass_and_load - quad mass"
    end if
    u = 0.0d0
    do i = 1, 4
        u(6*i-2:6*i) = axis
    end do
    if (abs(dot_product(u, matmul(m, u)) - rho * t**3 * area / 12.0d0) > &
        tol * rho * t**3 * area / 12.0d0) then
        rst = .false.
        print "(A)", "TEST FAILED: test_shell_mass_and_load - quad inertia"
    end if

    ! The consistent nodal loads sum to the total applied load
    f = quad%external_force_vector(q)
    ftot = 0.0d0
    do i = 1, 4
        ftot = ftot + f(6*i-5:6*i-3)
    end do
    if (maxval(abs(ftot - area * q)) > tol) then
        rst = .false.
        print "(A)", "TEST FAILED: test_shell_mass_and_load - quad load"
    end if
end function

! ------------------------------------------------------------------------------
function test_shell_cantilever_static() result(rst)
    ! A thin plate strip (nu = 0) clamped at x = 0 with a unit tip load
    ! should match Euler-Bernoulli beam theory: delta = P L^3 / (3 E I).
    logical :: rst
    integer(int32), parameter :: nx = 10
    integer(int32), parameter :: ny = 2
    real(real64), parameter :: len = 1.0d0
    real(real64), parameter :: b = 0.1d0
    real(real64), parameter :: t = 0.01d0
    real(real64), parameter :: e = 2.0d11
    real(real64), parameter :: rtol = 1.0d-2
    integer(int32) :: i, j, c, k1, k2, k3, k4
    real(real64) :: expected, delta
    type(material) :: mat
    type(node), allocatable, dimension(:) :: nodes
    type(rectangular_shell_element), dimension(nx * ny) :: quads
    type(triangular_shell_element), dimension(2 * nx * ny) :: tris

    rst = .true.
    mat = material(e, 0.0d0, 7.85d3)
    expected = len**3 / (3.0d0 * e * b * t**3 / 12.0d0)

    ! Quadrilateral mesh in the x-y plane, loaded along global z
    call build_strip(nx, ny, len, b, .false., nodes)
    c = 0
    do j = 1, ny
        do i = 1, nx
            c = c + 1
            k1 = (j - 1) * (nx + 1) + i
            k2 = k1 + 1
            k3 = k2 + nx + 1
            k4 = k1 + nx + 1
            quads(c) = rectangular_shell_element(mat, t, nodes(k1), &
                nodes(k2), nodes(k3), nodes(k4))
        end do
    end do
    delta = cantilever_tip_deflection(quads, nodes, nx, ny, 3)
    if (abs(delta - expected) > rtol * expected) then
        rst = .false.
        print "(A)", "TEST FAILED: test_shell_cantilever_static - quad"
        print "(A, ES12.5, A, ES12.5)", "Expected: ", expected, ", Found: ", delta
    end if

    ! Triangular mesh in the x-z plane, loaded along global y
    call build_strip(nx, ny, len, b, .true., nodes)
    c = 0
    do j = 1, ny
        do i = 1, nx
            k1 = (j - 1) * (nx + 1) + i
            k2 = k1 + 1
            k3 = k2 + nx + 1
            k4 = k1 + nx + 1
            c = c + 1
            tris(c) = triangular_shell_element(mat, t, nodes(k1), &
                nodes(k2), nodes(k3))
            c = c + 1
            tris(c) = triangular_shell_element(mat, t, nodes(k1), &
                nodes(k3), nodes(k4))
        end do
    end do
    delta = abs(cantilever_tip_deflection(tris, nodes, nx, ny, 2))
    if (abs(delta - expected) > rtol * expected) then
        rst = .false.
        print "(A)", "TEST FAILED: test_shell_cantilever_static - tri"
        print "(A, ES12.5, A, ES12.5)", "Expected: ", expected, ", Found: ", delta
    end if
end function

! ------------------------------------------------------------------------------
function test_shell_cantilever_modal() result(rst)
    ! The first natural frequency of a thin clamped plate strip (nu = 0)
    ! should match Euler-Bernoulli beam theory:
    ! omega = (1.8751)^2 sqrt(E I / (rho A L^4)).
    logical :: rst
    integer(int32), parameter :: nx = 10
    integer(int32), parameter :: ny = 2
    real(real64), parameter :: len = 1.0d0
    real(real64), parameter :: b = 0.1d0
    real(real64), parameter :: t = 0.01d0
    real(real64), parameter :: e = 2.0d11
    real(real64), parameter :: rho = 7.85d3
    real(real64), parameter :: rtol = 1.0d-2
    integer(int32) :: i, j, c, k1, k2, k3, k4
    real(real64) :: expected, omega
    type(material) :: mat
    type(node), allocatable, dimension(:) :: nodes
    type(rectangular_shell_element), dimension(nx * ny) :: quads
    type(triangular_shell_element), dimension(2 * nx * ny) :: tris

    rst = .true.
    mat = material(e, 0.0d0, rho)
    expected = 1.875104068711961d0**2 * &
        sqrt(e * b * t**3 / 12.0d0 / (rho * b * t * len**4))

    call build_strip(nx, ny, len, b, .false., nodes)
    c = 0
    do j = 1, ny
        do i = 1, nx
            k1 = (j - 1) * (nx + 1) + i
            k2 = k1 + 1
            k3 = k2 + nx + 1
            k4 = k1 + nx + 1
            c = c + 1
            quads(c) = rectangular_shell_element(mat, t, nodes(k1), &
                nodes(k2), nodes(k3), nodes(k4))
            tris(2*c-1) = triangular_shell_element(mat, t, nodes(k1), &
                nodes(k2), nodes(k3))
            tris(2*c) = triangular_shell_element(mat, t, nodes(k1), &
                nodes(k3), nodes(k4))
        end do
    end do

    omega = cantilever_first_frequency(quads, nodes, nx, ny)
    if (abs(omega - expected) > rtol * expected) then
        rst = .false.
        print "(A)", "TEST FAILED: test_shell_cantilever_modal - quad"
        print "(A, ES12.5, A, ES12.5)", "Expected: ", expected, ", Found: ", omega
    end if

    omega = cantilever_first_frequency(tris, nodes, nx, ny)
    if (abs(omega - expected) > rtol * expected) then
        rst = .false.
        print "(A)", "TEST FAILED: test_shell_cantilever_modal - tri"
        print "(A, ES12.5, A, ES12.5)", "Expected: ", expected, ", Found: ", omega
    end if
end function

! ------------------------------------------------------------------------------
end module

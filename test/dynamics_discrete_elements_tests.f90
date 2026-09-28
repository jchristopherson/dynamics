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
module dynamics_discrete_elements_tests
    use iso_fortran_env
    use dynamics
    use fortran_test_helper
    implicit none

contains
! ------------------------------------------------------------------------------
function test_spring_elements() result(rst)
    logical :: rst
    integer(int32) :: row, col
    real(real64), parameter :: tol = 1.0d-12
    real(real64), dimension(4) :: b2
    real(real64), dimension(6) :: b3
    real(real64), dimension(4,4) :: expected_k2
    real(real64), dimension(6,6) :: expected_k3
    real(real64), allocatable, dimension(:) :: elongation, force, s
    real(real64), allocatable, dimension(:,:) :: k, m, c, kdense
    type(csr_matrix) :: kcsr
    type(node), dimension(2) :: nodes_2d, nodes_3d, beam_nodes
    type(spring_element_2d), dimension(1) :: springs_2d
    type(spring_element_3d), dimension(1) :: springs_3d

    rst = .true.

    ! 2D spring with the axis defined by the nodes
    nodes_2d(1) = node(1, 2, 0.0d0, 0.0d0, 0.0d0)
    nodes_2d(2) = node(2, 2, 3.0d0, 4.0d0, 0.0d0)
    springs_2d(1) = spring_element_2d(10.0d0, nodes_2d(1), nodes_2d(2))
    b2 = [-0.6d0, -0.8d0, 0.6d0, 0.8d0]
    do col = 1, 4
        do row = 1, 4
            expected_k2(row,col) = 10.0d0 * b2(row) * b2(col)
        end do
    end do
    call assemble_static_system(4, springs_2d, nodes_2d, k)
    m = springs_2d(1)%mass_matrix()
    c = springs_2d(1)%damping_matrix()
    elongation = springs_2d(1)%strain([0.0d0, 0.0d0, 0.3d0, 0.4d0], [0.0d0])
    force = springs_2d(1)%stress([0.0d0, 0.0d0, 0.3d0, 0.4d0], [0.0d0])
    s = springs_2d(1)%get_node_natural_coordinates(2)
    if (maxval(abs(k - expected_k2)) > tol .or. &
        maxval(abs(m)) > 0.0d0 .or. size(m, 1) /= 4 .or. &
        maxval(abs(c)) > 0.0d0 .or. size(c, 1) /= 4 .or. &
        abs(elongation(1) - 0.5d0) > tol .or. &
        abs(force(1) - 5.0d0) > tol .or. &
        abs(s(1) - 1.0d0) > tol) then
        rst = .false.
        print "(A)", "TEST FAILED: test_spring_elements - 2D"
    end if

    ! 2D zero-length spring with a user-defined axis
    nodes_2d(2) = node(2, 2, 0.0d0, 0.0d0, 0.0d0)
    springs_2d(1) = spring_element_2d(10.0d0, nodes_2d(1), nodes_2d(2), &
        [0.0d0, 2.0d0])
    b2 = [0.0d0, -1.0d0, 0.0d0, 1.0d0]
    do col = 1, 4
        do row = 1, 4
            expected_k2(row,col) = 10.0d0 * b2(row) * b2(col)
        end do
    end do
    k = springs_2d(1)%stiffness_matrix()
    if (maxval(abs(k - expected_k2)) > tol) then
        rst = .false.
        print "(A)", "TEST FAILED: test_spring_elements - 2D user axis"
    end if

    ! 2D spring attached to nodes carrying a rotational DOF (e.g. beam nodes)
    beam_nodes(1) = node(1, 3, 0.0d0, 0.0d0, 0.0d0)
    beam_nodes(2) = node(2, 3, 2.0d0, 0.0d0, 0.0d0)
    springs_2d(1) = spring_element_2d(5.0d0, beam_nodes(1), beam_nodes(2))
    call assemble_static_system(6, springs_2d, beam_nodes, k)
    if (abs(k(1,1) - 5.0d0) > tol .or. abs(k(4,4) - 5.0d0) > tol .or. &
        abs(k(1,4) + 5.0d0) > tol .or. abs(k(4,1) + 5.0d0) > tol .or. &
        abs(sum(abs(k)) - 20.0d0) > tol) then
        rst = .false.
        print "(A)", "TEST FAILED: test_spring_elements - 2D beam nodes"
    end if

    ! 3D spring with the axis defined by the nodes
    nodes_3d(1) = node(1, 3, 0.0d0, 0.0d0, 0.0d0)
    nodes_3d(2) = node(2, 3, 2.0d0, 3.0d0, 6.0d0)
    springs_3d(1) = spring_element_3d(14.0d0, nodes_3d(1), nodes_3d(2))
    b3 = [-2.0d0, -3.0d0, -6.0d0, 2.0d0, 3.0d0, 6.0d0] / 7.0d0
    do col = 1, 6
        do row = 1, 6
            expected_k3(row,col) = 14.0d0 * b3(row) * b3(col)
        end do
    end do
    call assemble_static_system(6, springs_3d, nodes_3d, k)
    call assemble_static_system(6, springs_3d, nodes_3d, kcsr)
    allocate(kdense(6,6))
    kdense = kcsr
    elongation = springs_3d(1)%strain([0.0d0, 0.0d0, 0.0d0, &
        2.0d0, 3.0d0, 6.0d0], [0.0d0])
    if (maxval(abs(k - expected_k3)) > tol .or. &
        maxval(abs(kdense - expected_k3)) > tol .or. &
        abs(elongation(1) - 7.0d0) > tol) then
        rst = .false.
        print "(A)", "TEST FAILED: test_spring_elements - 3D"
    end if

    ! 3D zero-length spring with a user-defined axis
    nodes_3d(2) = node(2, 3, 0.0d0, 0.0d0, 0.0d0)
    springs_3d(1) = spring_element_3d(14.0d0, nodes_3d(1), nodes_3d(2), &
        [0.0d0, 0.0d0, 3.0d0])
    b3 = [0.0d0, 0.0d0, -1.0d0, 0.0d0, 0.0d0, 1.0d0]
    do col = 1, 6
        do row = 1, 6
            expected_k3(row,col) = 14.0d0 * b3(row) * b3(col)
        end do
    end do
    k = springs_3d(1)%stiffness_matrix()
    if (maxval(abs(k - expected_k3)) > tol) then
        rst = .false.
        print "(A)", "TEST FAILED: test_spring_elements - 3D user axis"
    end if
end function

! ------------------------------------------------------------------------------
function test_damper_elements() result(rst)
    logical :: rst
    integer(int32) :: row, col
    real(real64), parameter :: tol = 1.0d-12
    real(real64), dimension(4) :: b2
    real(real64), dimension(6) :: b3
    real(real64), dimension(4,4) :: expected_c2
    real(real64), dimension(6,6) :: expected_c3
    real(real64), allocatable, dimension(:) :: force
    real(real64), allocatable, dimension(:,:) :: c, k, m, cdense
    type(csr_matrix) :: ccsr
    type(node), dimension(2) :: nodes_2d, nodes_3d
    type(damper_element_2d), dimension(1) :: dampers_2d
    type(damper_element_3d), dimension(1) :: dampers_3d

    rst = .true.

    ! 2D damper
    nodes_2d(1) = node(1, 2, 0.0d0, 0.0d0, 0.0d0)
    nodes_2d(2) = node(2, 2, 0.0d0, 2.0d0, 0.0d0)
    dampers_2d(1) = damper_element_2d(4.0d0, nodes_2d(1), nodes_2d(2))
    b2 = [0.0d0, -1.0d0, 0.0d0, 1.0d0]
    do col = 1, 4
        do row = 1, 4
            expected_c2(row,col) = 4.0d0 * b2(row) * b2(col)
        end do
    end do
    call assemble_damping_matrix(4, dampers_2d, nodes_2d, c)
    k = dampers_2d(1)%stiffness_matrix()
    m = dampers_2d(1)%mass_matrix()
    force = dampers_2d(1)%stress([0.0d0, 0.0d0, 0.0d0, 0.5d0], [0.0d0])
    if (maxval(abs(c - expected_c2)) > tol .or. &
        maxval(abs(k)) > 0.0d0 .or. size(k, 1) /= 4 .or. &
        maxval(abs(m)) > 0.0d0 .or. size(m, 1) /= 4 .or. &
        abs(force(1) - 2.0d0) > tol) then
        rst = .false.
        print "(A)", "TEST FAILED: test_damper_elements - 2D"
    end if

    ! 3D damper
    nodes_3d(1) = node(1, 3, 0.0d0, 0.0d0, 0.0d0)
    nodes_3d(2) = node(2, 3, 1.0d0, 2.0d0, 2.0d0)
    dampers_3d(1) = damper_element_3d(9.0d0, nodes_3d(1), nodes_3d(2))
    b3 = [-1.0d0, -2.0d0, -2.0d0, 1.0d0, 2.0d0, 2.0d0] / 3.0d0
    do col = 1, 6
        do row = 1, 6
            expected_c3(row,col) = 9.0d0 * b3(row) * b3(col)
        end do
    end do
    call assemble_damping_matrix(6, dampers_3d, nodes_3d, c)
    call assemble_damping_matrix(6, dampers_3d, nodes_3d, ccsr)
    allocate(cdense(6,6))
    cdense = ccsr
    k = dampers_3d(1)%stiffness_matrix()
    if (maxval(abs(c - expected_c3)) > tol .or. &
        maxval(abs(cdense - expected_c3)) > tol .or. &
        maxval(abs(k)) > 0.0d0 .or. size(k, 1) /= 6) then
        rst = .false.
        print "(A)", "TEST FAILED: test_damper_elements - 3D"
    end if
end function

! ------------------------------------------------------------------------------
function test_mass_elements() result(rst)
    logical :: rst
    integer(int32) :: i
    real(real64), parameter :: tol = 1.0d-12
    real(real64), dimension(6,6) :: expected_m
    real(real64), allocatable, dimension(:) :: s
    real(real64), allocatable, dimension(:,:) :: m, k, c
    type(node), dimension(1) :: nodes_2d, nodes_3d
    type(mass_element_2d), dimension(1) :: masses_2d
    type(mass_element_3d), dimension(1) :: masses_3d

    rst = .true.

    ! 2D point mass
    nodes_2d(1) = node(1, 2, 1.0d0, 2.0d0, 0.0d0)
    masses_2d(1) = mass_element_2d(3.0d0, nodes_2d(1))
    call assemble_dynamic_system(2, masses_2d, nodes_2d, m, k)
    c = masses_2d(1)%damping_matrix()
    s = masses_2d(1)%get_node_natural_coordinates(1)
    if (abs(m(1,1) - 3.0d0) > tol .or. abs(m(2,2) - 3.0d0) > tol .or. &
        abs(m(1,2)) > tol .or. abs(m(2,1)) > tol .or. &
        maxval(abs(k)) > 0.0d0 .or. maxval(abs(c)) > 0.0d0 .or. &
        size(c, 1) /= 2 .or. abs(s(1)) > tol) then
        rst = .false.
        print "(A)", "TEST FAILED: test_mass_elements - 2D"
    end if

    ! 3D point mass attached to a node carrying rotational DOFs
    nodes_3d(1) = node(1, 6, 0.0d0, 0.0d0, 0.0d0)
    masses_3d(1) = mass_element_3d(2.0d0, nodes_3d(1))
    expected_m = 0.0d0
    do i = 1, 3
        expected_m(i,i) = 2.0d0
    end do
    call assemble_dynamic_system(6, masses_3d, nodes_3d, m, k)
    if (maxval(abs(m - expected_m)) > tol .or. maxval(abs(k)) > 0.0d0) then
        rst = .false.
        print "(A)", "TEST FAILED: test_mass_elements - 3D"
    end if
end function

! ------------------------------------------------------------------------------
function test_discrete_element_system() result(rst)
    ! Single DOF spring-mass-damper system constructed from discrete elements
    logical :: rst
    real(real64), parameter :: tol = 1.0d-10
    real(real64), parameter :: stiffness = 100.0d0
    real(real64), parameter :: damping = 2.0d0
    real(real64), parameter :: mass = 4.0d0
    real(real64), parameter :: load = 10.0d0
    integer(int32), dimension(3) :: constraints
    real(real64), dimension(4) :: f
    real(real64), allocatable, dimension(:) :: fr, u
    real(real64), allocatable, dimension(:,:) :: k, c, m, kunused, kr, cr, mr
    type(node), dimension(2) :: nodes
    type(spring_element_2d), dimension(1) :: springs
    type(damper_element_2d), dimension(1) :: dampers
    type(mass_element_2d), dimension(1) :: masses

    rst = .true.
    nodes(1) = node(1, 2, 0.0d0, 0.0d0, 0.0d0)
    nodes(2) = node(2, 2, 1.0d0, 0.0d0, 0.0d0)
    springs(1) = spring_element_2d(stiffness, nodes(1), nodes(2))
    dampers(1) = damper_element_2d(damping, nodes(1), nodes(2))
    masses(1) = mass_element_2d(mass, nodes(2))

    call assemble_static_system(4, springs, nodes, k)
    call assemble_damping_matrix(4, dampers, nodes, c)
    call assemble_dynamic_system(4, masses, nodes, m, kunused)

    ! Fix node 1 entirely and the y-translation of node 2
    constraints = [1, 2, 4]
    kr = apply_boundary_conditions(constraints, k)
    cr = apply_boundary_conditions(constraints, c)
    mr = apply_boundary_conditions(constraints, m)
    f = [0.0d0, 0.0d0, load, 0.0d0]
    fr = apply_boundary_conditions(constraints, f)
    u = solve_static_system(kr, fr)

    if (size(kr, 1) /= 1 .or. &
        abs(kr(1,1) - stiffness) > tol .or. &
        abs(cr(1,1) - damping) > tol .or. &
        abs(mr(1,1) - mass) > tol .or. &
        abs(u(1) - load / stiffness) > tol .or. &
        abs(sqrt(kr(1,1) / mr(1,1)) - 5.0d0) > tol) then
        rst = .false.
        print "(A)", "TEST FAILED: test_discrete_element_system"
    end if
end function

! ------------------------------------------------------------------------------
function test_assemble_discrete_system() result(rst)
    ! Two mass-spring-damper pairs in series
    logical :: rst
    real(real64), parameter :: tol = 1.0d-12
    real(real64), allocatable, dimension(:,:) :: m, c, k, mref, cref, kref, &
        kunused, mdense, cdense, kdense
    type(csr_matrix) :: mcsr, ccsr, kcsr
    type(node), dimension(3) :: nodes
    type(spring_element_2d), dimension(2) :: springs
    type(damper_element_2d), dimension(2) :: dampers
    type(mass_element_2d), dimension(2) :: masses
    type(damper_element_2d), dimension(0) :: no_dampers

    rst = .true.
    nodes(1) = node(1, 2, 0.0d0, 0.0d0, 0.0d0)
    nodes(2) = node(2, 2, 0.0d0, 1.0d0, 0.0d0)
    nodes(3) = node(3, 2, 0.0d0, 2.0d0, 0.0d0)
    springs(1) = spring_element_2d(100.0d0, nodes(1), nodes(2))
    springs(2) = spring_element_2d(200.0d0, nodes(2), nodes(3))
    dampers(1) = damper_element_2d(3.0d0, nodes(1), nodes(2))
    dampers(2) = damper_element_2d(5.0d0, nodes(2), nodes(3))
    masses(1) = mass_element_2d(1.5d0, nodes(2))
    masses(2) = mass_element_2d(2.5d0, nodes(3))

    ! Reference matrices from the individual assembly routines
    call assemble_dynamic_system(6, masses, nodes, mref, kunused)
    call assemble_static_system(6, springs, nodes, kref)
    call assemble_damping_matrix(6, dampers, nodes, cref)

    call assemble_discrete_system(6, masses, dampers, springs, nodes, m, c, k)
    if (maxval(abs(m - mref)) > tol .or. maxval(abs(c - cref)) > tol .or. &
        maxval(abs(k - kref)) > tol) then
        rst = .false.
        print "(A)", "TEST FAILED: test_assemble_discrete_system - dense"
    end if

    call assemble_discrete_system(6, masses, dampers, springs, nodes, &
        mcsr, ccsr, kcsr)
    allocate(mdense(6,6), cdense(6,6), kdense(6,6))
    mdense = mcsr
    cdense = ccsr
    kdense = kcsr
    if (maxval(abs(mdense - mref)) > tol .or. &
        maxval(abs(cdense - cref)) > tol .or. &
        maxval(abs(kdense - kref)) > tol) then
        rst = .false.
        print "(A)", "TEST FAILED: test_assemble_discrete_system - CSR"
    end if

    call assemble_discrete_system(6, masses, no_dampers, springs, nodes, &
        m, c, k)
    if (maxval(abs(m - mref)) > tol .or. maxval(abs(c)) > 0.0d0 .or. &
        maxval(abs(k - kref)) > tol) then
        rst = .false.
        print "(A)", "TEST FAILED: test_assemble_discrete_system - no dampers"
    end if
end function

! ------------------------------------------------------------------------------
end module

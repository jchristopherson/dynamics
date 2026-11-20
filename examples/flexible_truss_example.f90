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
program example
    use iso_fortran_env, only : real64
    use dynamics, only : apply_boundary_conditions, resolve_boundary_conditions, nodal_condition, &
        DOF_UX, DOF_UY, DOF_UZ, DOF_RX, DOF_RY, DOF_RZ, &
        euler_bernoulli_beam_element, fe_node, finite_element_mesh, &
        flexible_material, linear_elastic_isotropic_material
    use fplot_core, only : plot_2d, plot_data_2d
    use linalg, only : csr_matrix, size
    implicit none

    ! Structure Geometry (same as structural_example):
    !
    !         2
    !        /|\
    !       / | \
    !      /  |  \
    !     1---4---3
    !
    ! All nodes lie in the XY plane, but the spatial beams give each node six
    ! DOFs: three translations and three rotations. This example assembles the
    ! reference mesh and applies the supports without loads or solves.

    ! Mesh Refinement
    integer, parameter :: nel = 50
        !! Number of Euler-Bernoulli elements per leg.
    integer, parameter :: nlegs = 5
    integer, parameter :: ncorners = 4
    integer, parameter :: nnodes = ncorners + nlegs * (nel - 1)
    integer, parameter :: nelements = nlegs * nel

    ! Corner Coordinates (m) and Leg Connectivity
    real(real64), parameter, dimension(ncorners) :: corner_x = [0.0d0, 0.5d0, 1.0d0, 0.5d0]
    real(real64), parameter, dimension(ncorners) :: corner_y = [0.0d0, 0.2d0, 0.0d0, 0.0d0]
    integer, parameter, dimension(nlegs) :: leg_start = [1, 2, 3, 4, 2]
    integer, parameter, dimension(nlegs) :: leg_end = [2, 3, 4, 1, 4]
    real(real64), parameter, dimension(3) :: reference_y = [0.0d0, 0.0d0, 1.0d0]
        !! Normal to the truss plane, so it cannot be parallel to any beam axis.

    ! Solid Circular Section and Material (SI units)
    real(real64), parameter :: pi = acos(-1.0d0)
    real(real64), parameter :: diameter = 5.0d-2
    real(real64), parameter :: area = pi * diameter**2 / 4.0d0
    real(real64), parameter :: moi = pi * diameter**4 / 64.0d0
    real(real64), parameter :: torsion_constant = 2.0d0 * moi
    real(real64), parameter :: modulus = 70.0d9
    real(real64), parameter :: poissons_ratio = 0.33d0
    real(real64), parameter :: density = 2.7d3

    type(linear_elastic_isotropic_material), target :: material
    class(flexible_material), pointer :: mat
    type(finite_element_mesh) :: mesh
    type(csr_matrix) :: mass
    type(csr_matrix) :: stiffness
    integer, allocatable, dimension(:) :: bc
    real(real64), allocatable, dimension(:) :: q
    integer, dimension(nel + 1, nlegs) :: leg_nodes
    integer :: ndof
    integer :: i
    integer :: leg
    integer :: el
    integer :: node_count
    integer :: elem_count
    real(real64) :: t
    real(real64), dimension(3) :: start
    real(real64), dimension(3) :: finish

    ! The elements share a pointer to this material, which outlives the mesh.
    material = linear_elastic_isotropic_material(modulus, poissons_ratio, density)
    mat => material

    mesh = finite_element_mesh(nnodes, nelements)
    do i = 1, ncorners
        mesh%nodes(i) = fe_node(i, [corner_x(i), corner_y(i), 0.0d0])
    end do

    ! Reuse nodes 1-4 at the junctions so meeting legs share DOFs. Number only
    ! the interior nodes as new nodes, interpolating their coordinates along
    ! each leg. Keep the ordered node indices for both assembly and plotting.
    node_count = ncorners
    elem_count = 0
    do leg = 1, nlegs
        leg_nodes(1, leg) = leg_start(leg)
        leg_nodes(nel + 1, leg) = leg_end(leg)
        start = mesh%nodes(leg_start(leg))%coordinates
        finish = mesh%nodes(leg_end(leg))%coordinates

        do i = 1, nel - 1
            t = real(i, real64) / real(nel, real64)
            node_count = node_count + 1
            mesh%nodes(node_count) = fe_node(node_count, start + t * (finish - start))
            leg_nodes(i + 1, leg) = node_count
        end do

        do el = 1, nel
            elem_count = elem_count + 1
            call mesh%assign_element(elem_count, &
                euler_bernoulli_beam_element( &
                    mat, &
                    area, &
                    torsion_constant, &
                    moi, &
                    moi, &
                    leg_nodes(el, leg), &
                    leg_nodes(el + 1, leg), &
                    reference_y &
                ) &
            )
        end do
    end do

    ! Fix all six DOFs at node 1. At node 3, restrain only the vertical
    ! translation, as in structural_example.
    call resolve_boundary_conditions(mesh, [ &
        nodal_condition(1, DOF_UX), nodal_condition(1, DOF_UY), nodal_condition(1, DOF_UZ), &
        nodal_condition(1, DOF_RX), nodal_condition(1, DOF_RY), nodal_condition(1, DOF_RZ), &
        nodal_condition(3, DOF_UY)], bc, ndof)

    ! The zero state selects the undeformed tangent; the CSR outputs select
    ! direct sparse assembly. Use the same support equations for both matrices.
    allocate(q(ndof), source=0.0d0)
    call mesh%assemble_mass(mass)
    call mesh%assemble_stiffness(q, stiffness)

    mass = apply_boundary_conditions(bc, mass)
    stiffness = apply_boundary_conditions(bc, stiffness)

    print "(A,I0,A,I0,A,I0,A)", "Mesh: ", size(mesh%nodes), " nodes, ", size(mesh%elements), &
        " elements, ", size(q), " DOFs"
    print "(A,I0)", "Free DOFs after applying supports: ", size(stiffness, 1)
    call plot_mesh()

contains
! ------------------------------------------------------------------------------
    subroutine plot_mesh()
        !! Draw the reference centerlines with a marker at each mesh node.
        type(plot_2d) :: plt
        type(plot_data_2d), dimension(nlegs) :: data
        real(real64), dimension(nel + 1) :: x
        real(real64), dimension(nel + 1) :: y
        integer :: leg
        integer :: i
        integer :: nid

        call plt%initialize()
        call plt%set_title("Undeformed Euler-Bernoulli Truss Mesh")
        call plt%set_x_axis_title("X (m)")
        call plt%set_y_axis_title("Y (m)")
        call plt%show_legend(.false.)

        ! Keep a separate data series for each leg so FPLOT retains all five
        ! centerlines without drawing lines between unrelated endpoints.
        do leg = 1, nlegs
            do i = 1, nel + 1
                nid = leg_nodes(i, leg)
                x(i) = mesh%nodes(nid)%coordinates(1)
                y(i) = mesh%nodes(nid)%coordinates(2)
            end do
            call data(leg)%define_data(x, y)
            call data(leg)%set_draw_line(.true.)
            call data(leg)%set_draw_markers(.true.)
            call data(leg)%set_marker_scaling(0.35)
            call plt%push(data(leg))
        end do

        call plt%draw()
    end subroutine
end program
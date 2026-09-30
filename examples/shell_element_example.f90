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
    !! This example considers a rectangular shell with fixed supports at
    !! opposing edges.  The analysis will compute the first 6 mode shapes and 
    !! plot the results.
    use iso_fortran_env
    use dynamics
    use fplot_core
    use linalg
    implicit none

    ! Material Properties
    real(real64), parameter :: modulus = 70.0d9
    real(real64), parameter :: density = 2.7d3
    real(real64), parameter :: poissons_ratio = 0.33d0

    ! Geometry
    real(real64), parameter :: thickness = 1.0d-2
    real(real64), parameter :: width = 0.75d0
    real(real64), parameter :: length = 1.25d0
    integer(int32), parameter :: ndiv = 50
    integer(int32), parameter :: dof_per_node = 6
    integer(int32), parameter :: number_of_nodes = ndiv * ndiv
    integer(int32), parameter :: global_dof_count = dof_per_node * number_of_nodes
    integer(int32), parameter :: number_of_elements = (ndiv - 1)**2
    integer(int32), parameter :: number_of_boundary_conditions = 2 * dof_per_node * ndiv

    ! Analysis Properties
    integer(int32), parameter :: mode_count = 6
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    integer(int32) :: i, j, kk, n1, n2, n3, n4, n5, n6
    integer(int32), allocatable, dimension(:) :: bc
    real(real64), allocatable, dimension(:) :: xc, yc, freqs
    real(real64), allocatable, dimension(:,:) :: shapes, shapesbc, u
    real(real64), allocatable, dimension(:,:,:) :: xy
    type(material) :: mat
    type(node), allocatable, dimension(:) :: nodes
    type(rectangular_shell_element), allocatable, dimension(:) :: elements
    type(csr_matrix) :: K, Kbc, M, Mbc
    

! ------------------------------------------------------------------------------
! Create the mesh
    xc = linspace(0.0d0, width, ndiv)
    yc = linspace(0.0d0, length, ndiv)
    xy = meshgrid(xc, yc)               ! The x-y coordinates of each node

    ! Number the nodes.  The nodes will be numbered along the x-axis, column
    ! by column.
    allocate(nodes(number_of_nodes))
    kk = 0
    do j = 1, ndiv
        do i = 1, ndiv
            kk = kk + 1
            nodes(kk) = node(kk, dof_per_node, xy(i,j,1), xy(i,j,2), 0.0d0)
        end do
    end do

    ! Define the material
    mat = material(modulus, poissons_ratio, density)

    ! Build the element list
    allocate(elements(number_of_elements))
    kk = 0
    n1 = 1
    n2 = 2
    do j = 1, ndiv - 1
        do i = 1, ndiv - 1
            kk = kk + 1
            n3 = n2 + ndiv
            n4 = n3 - 1
            elements(kk) = rectangular_shell_element( &
                mat, &
                thickness, &
                nodes(n1), &
                nodes(n2), &
                nodes(n3), &
                nodes(n4) &
            )
            if (i == ndiv - 1) then
                n1 = n3 + 1 - ndiv
            else
                n1 = n2
            end if
            n2 = n1 + 1
        end do
    end do

! ------------------------------------------------------------------------------
! Build the FE Model

    ! Assemble the global matrices
    call assemble_dynamic_system(global_dof_count, elements, nodes, M, K)

    ! Define the boundary conditions
    allocate(bc(number_of_boundary_conditions))
    ! First Edge
    do i = 1, number_of_boundary_conditions / 2
        bc(i) = i
    end do
    kk = dof_per_node * (ndiv * (ndiv - 1) + 1) - dof_per_node + 1
    do i = number_of_boundary_conditions / 2 + 1, number_of_boundary_conditions
        bc(i) = kk
        kk = kk + 1
    end do

    ! Apply the boundary conditions
    Mbc = apply_boundary_conditions(bc, M)
    Kbc = apply_boundary_conditions(bc, K)

! ------------------------------------------------------------------------------
! Analysis

    ! Compute the modal solution
    call modal_response(Mbc, Kbc, mode_count, freqs, shapesbc)

    ! Conver the frequency units to Hz, from rad/s
    freqs = freqs / (2.0d0 * pi)

    ! Reconstruct each mode shape, including the boundary conditions
    allocate(shapes(size(K, 1), mode_count))
    do i = 1, mode_count
        shapes(:,i) = restore_constrained_values(bc, shapesbc(:,i))
    end do

    ! Write out each resonant frequency
    print "(A)", "MODAL RESPONSE:"
    do i = 1, size(freqs)
        print "(A, I0, A, F0.3, A)", "Mode ", i, ": ", freqs(i), " Hz"
    end do

! ------------------------------------------------------------------------------
! Plots
    allocate(u(ndiv, ndiv))
    do i = 1, size(freqs)
        call extract_nodal_displacements(shapes(:,i), u)
        call plot_mode_shape(i, freqs(i), xy(:,:,1), xy(:,:,2), u)
    end do

contains
! ------------------------------------------------------------------------------
    subroutine plot_mode_shape(mode, freq, x, y, shape)
        !! Plots the specified mode shape.
        integer(int32), intent(in) :: mode
            !! The mode number.
        real(real64), intent(in) :: freq
            !! The frequency, in Hz.
        real(real64), intent(in), dimension(:,:) :: x
            !! The mesh x-coordinates.
        real(real64), intent(in), dimension(:,:) :: y
            !! The mesh y-coordinates.
        real(real64), intent(in), dimension(:,:) :: shape
            !! The mode shape corresponding the x-y points.

        ! Local Variables
        type(surface_plot) :: plt
        type(rainbow_colormap) :: map
        character(len = 40) :: title

        ! Define the title
        write (title, "(A, I0, A, F0.3, A)") "Mode ", mode, ": ", freq, " Hz"

        ! Plot the results
        call plt%initialize()
        call plt%set_elevation(40.0d0)
        call plt%set_azimuth(50.0d0)
        call plt%set_colormap(map)
        call plt%set_title(trim(title))
        call plt%set_x_axis_title("x")
        call plt%set_y_axis_title("y")
        call plt%push(x, y, shape)
        call plt%draw()
    end subroutine

! ------------------------------------------------------------------------------
    subroutine extract_nodal_displacements(shape, du)
        !! Extracts the x, y, and z translations of each node from a mode
        !! shape vector and arranges them on the meshgrid layout such that
        !! entry (i,j) corresponds to the point (xy(i,j,1), xy(i,j,2)).
        real(real64), intent(in), dimension(:) :: shape
            !! The global mode shape vector, including constrained DOF.
        real(real64), intent(out), dimension(:,:) :: du
            !! The ndiv-by-ndiv matrix of displacements.

        ! Local Variables
        integer(int32) :: ii, jj, dof

        ! Process
        do jj = 1, ndiv
            do ii = 1, ndiv
                ! Matches the node numbering used when building the mesh
                dof = dof_per_node * ((jj - 1) * ndiv + ii - 1)
                du(ii,jj) = shape(dof + 3)
            end do
        end do
    end subroutine

! ------------------------------------------------------------------------------
end program
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
    use iso_fortran_env
    use dynamics
    use linalg
    implicit none

    ! Structure Geometry:
    !
    ! The structure is a simple, 2D truss structure shaped as follows.
    !
    !       2
    !      /|\
    !     / | \
    !    /  |  \
    !   /   |   \
    ! 1 ---------- 3
    !       4
    !
    ! The truss is pinned from both translational DOF at it's lower left corner,
    ! and supported vertically at it's lower right corner.  A static analysis
    ! is performed by applying a vertical force at the apex of the structure
    ! acting downwards.  Each leg of the structure (1-2, 2-3, 3-4, 4-1, and 2-4)
    ! is subdivided into NEL 2D beam elements, giving a total of 5 * NEL
    ! elements.  Because of the small size of this problem, dense matrices will
    ! be utilized.
    !
    ! Each node has 3 degrees-of-freedom: the two in-plane translations along 
    ! with one rotation.  The 4 corner nodes retain their original numbering
    ! (1 through 4); the interior nodes introduced by subdividing each leg are
    ! numbered afterwards, leg by leg.

    ! Constants
    character, parameter :: tab = achar(9)
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Mesh Refinement
    integer(int32), parameter :: nel = 50
        !! The number of beam elements used to represent each leg of the truss.
    integer(int32), parameter :: nlegs = 5
        !! The number of legs in the truss.
    integer(int32), parameter :: ncorners = 4
        !! The number of corner nodes in the truss.
    integer(int32), parameter :: dof_per_node = 3
    integer(int32), parameter :: nnodes = ncorners + nlegs * (nel - 1)
    integer(int32), parameter :: nelements = nlegs * nel
    integer(int32), parameter :: gdof = nnodes * dof_per_node
    integer(int32), parameter :: n_plot_modes = 6
        !! The number of lowest-frequency mode shapes to plot; the finer mesh
        !! yields many more modes than are of interest to visualize.

    ! Corner Node Coordinates & Leg Connectivity (1-2, 2-3, 3-4, 4-1, 2-4)
    real(real64), parameter :: corner_x(ncorners) = [0.0d0, 0.5d0, 1.0d0, 0.5d0]
    real(real64), parameter :: corner_y(ncorners) = [0.0d0, 0.2d0, 0.0d0, 0.0d0]
    integer(int32), parameter :: leg_start(nlegs) = [1, 2, 3, 4, 2]
    integer(int32), parameter :: leg_end(nlegs) = [2, 3, 4, 1, 4]

    ! Applied Forces
    real(real64), parameter :: applied_force = 1.0d4

    ! Beam Cross-Sectional Properties
    real(real64), parameter :: diameter = 5.0d-2
    real(real64), parameter :: area = 0.25d0 * pi * diameter**2
    real(real64), parameter :: moi = pi * diameter**4 / 64.0d0

    ! Material Properties
    real(real64), parameter :: modulus = 70.0d9
    real(real64), parameter :: poissons_ratio = 0.33d0
    real(real64), parameter :: density = 2.7d3

    ! FRF Parameters
    integer(int32), parameter :: nfreq = 1000
    real(real64), parameter :: alpha = 1.0d-3
    real(real64), parameter :: beta = 2.0d-6
    real(real64), parameter :: minfreq = 2.0d0 * pi * 1.0d1
    real(real64), parameter :: maxfreq = 2.0d0 * pi * 1.0d3

    ! Variables
    type(material) :: mat
    type(node) :: nodes(nnodes)
    type(beam_element_2d) :: beams(nelements)
    integer(int32) :: i, el, leg, node_count, elem_count, bc(4)
    integer(int32) :: leg_nodes(nlegs, nel + 1)
    real(real64) :: t, x0, y0, x1v, y1v
    real(real64), allocatable, dimension(:) :: F, Fbc, ubc, u, Fr, freqs
    real(real64), allocatable, dimension(:,:) :: shapesbc, shapes
    type(csr_matrix) :: K, Kbc, M, Mbc
    type(frf) :: frsp
    procedure(modal_excite), pointer :: excitefcn
    character(len = 32) :: str

    ! Create the material
    mat = material(modulus, poissons_ratio, density)

    ! Create the corner nodes
    do i = 1, ncorners
        nodes(i) = node(i, dof_per_node, corner_x(i), corner_y(i), 0.0d0)
    end do

    ! Subdivide each leg into NEL elements, creating the interior nodes and
    ! elements along the way.  LEG_NODES records the global node index of
    ! every node along each leg, in order, for later use when plotting.
    node_count = ncorners
    elem_count = 0
    do leg = 1, nlegs
        leg_nodes(leg, 1) = leg_start(leg)
        leg_nodes(leg, nel + 1) = leg_end(leg)
        x0 = corner_x(leg_start(leg))
        y0 = corner_y(leg_start(leg))
        x1v = corner_x(leg_end(leg))
        y1v = corner_y(leg_end(leg))

        do el = 1, nel - 1
            t = real(el, real64) / real(nel, real64)
            node_count = node_count + 1
            nodes(node_count) = node(node_count, dof_per_node, &
                x0 + t * (x1v - x0), y0 + t * (y1v - y0), 0.0d0)
            leg_nodes(leg, el + 1) = node_count
        end do

        do el = 1, nel
            elem_count = elem_count + 1
            beams(elem_count) = beam_element_2d(mat, area, moi, &
                nodes(leg_nodes(leg, el)), nodes(leg_nodes(leg, el + 1)))
        end do
    end do

    ! Assemble the stiffness and mass matrices
    call assemble_dynamic_system(gdof, beams, nodes, M, K)

    ! The external forcing vector can be easily constructed noting the force
    ! acts on node 2 alone, and acts in the negative y direction.
    allocate(F(size(K, 1)), source = 0.0d0)
    F(5) = -applied_force

    ! With both the external force vector and the stiffness matrix assembled,
    ! we can move on to applying boundary conditions.  The boundary conditions
    ! are applied to node 1 (all 3 DOF), and node 3 (vertical DOF).  We 
    ! construct an array with the appropriate global DOF identifiers.
    bc = [1, 2, 3, 8]

    ! Now, we can apply the boundary conditions to both the stiffness matrix
    ! and the external forcing vector
    Fbc = apply_boundary_conditions(bc, F)
    Kbc = apply_boundary_conditions(bc, K)
    Mbc = apply_boundary_conditions(bc, M)

    ! Solve the system
    ubc = solve_static_system(Kbc, Fbc)

    ! Account for the boundary conditions.  This simply inserts zeros where
    ! the boundary conditions have been applied such that it's now easy to
    ! determine the reaction loads at each constraint.
    u = restore_constrained_values(bc, ubc) * 1.0d3 ! applying scaling factor for ease of display

! --------------------
    ! Write out the displacement results for each node
    print "(A)", "NODAL DISPLACEMENTS (x 1000)"
    print "(A)", "Node" // tab // "   X" // tab // tab // "  Y" // tab  // tab // "   Theta"
    print "(A, EN0.3, A, EN0.3, A, EN0.3)", "1" // tab, u(1), tab, u(2), tab, u(3)
    print "(A, EN0.3, A, EN0.3, A, EN0.3)", "2" // tab, u(4), tab, u(5), tab, u(6)
    print "(A, EN0.3, A, EN0.3, A, EN0.3)", "3" // tab, u(7), tab, u(8), tab, u(9)
    print "(A, EN0.3, A, EN0.3, A, EN0.3)", "4" // tab, u(10), tab, u(11), tab, u(12)

    ! Compute the reaction loads
    u = 1.0d-3 * u  ! remove the scaling factor from above
    Fr = matmul(K, u)

    ! Write out the reaction loads at nodes 1 & 3
    print "(A)", new_line('a') // "REACTION LOADS"
    print "(A)", "Node" // tab // "  FX" // tab // tab // "  FY" // tab  // tab // "   MZ"
    print "(A, EN0.3, A, EN0.3, A, EN0.3)", "1" // tab, Fr(1), tab, Fr(2), tab, Fr(3)
    print "(A, EN0.3, A, EN0.3, A, EN0.3)", "3" // tab, Fr(7), tab, Fr(8), tab, Fr(9)
    print "(A, EN0.3, A, EN0.3, A, EN0.3)", "SUM" // tab, Fr(1) + Fr(7), tab, Fr(2) + Fr(8), tab, Fr(3) + Fr(9)

    ! Plot the deformed shape
    call plot_deflection(nodes, u, leg_nodes, 1.0d2, "Static Deformation")

! ******************************************************************************
! EIGEN ANALYSIS
! ------------------------------------------------------------------------------
    ! Compute the resonant modes of the system
    call modal_response(Mbc, Kbc, n_plot_modes, freqs, shapesbc)

    ! Convert the frequency units to units of Hz, from rad/s
    freqs = freqs / (2.0d0 * pi)

    ! Reconstruct each mode shape, including the boundary conditions
    allocate(shapes(size(u), size(shapesbc, 2)))
    do i = 1, size(shapesbc, 2)
        shapes(:,i) = restore_constrained_values(bc, shapesbc(:,i))
    end do

    ! Write out each mode shape
    print "(A)", new_line('a') // "MODAL RESPONSE"
    do i = 1, min(n_plot_modes, size(freqs))
        print "(A, I0, A, F0.3, A)", "Mode ", i, ": ", freqs(i), " Hz"
    end do

    ! Plot each of the lowest-frequency mode shapes
    do i = 1, min(n_plot_modes, size(freqs))
        write(str, "(A, I0, A, F0.3, A)") "Mode ", i, ": ", freqs(i), " Hz"
        call plot_deflection(nodes, shapes(:,i), leg_nodes, 1.0d-1, trim(str))
    end do

! ******************************************************************************
! HARMONIC ANALYSIS
! ------------------------------------------------------------------------------
    excitefcn => modal_frf_forcing_term
    frsp = frequency_response(Mbc, Kbc, alpha, beta, n_plot_modes, nfreq, &
        minfreq, maxfreq, excitefcn)
    call plot_frf(frsp, 5)

contains
! ------------------------------------------------------------------------------
    subroutine modal_frf_forcing_term(freq_, f_, args_)
        !! The forcing function.
        real(real64), intent(in) :: freq_
            !! The frequency value.
        complex(real64), intent(out), dimension(:) :: f_
            !! The output forcing vector.
        class(*), intent(inout), optional :: args_

        complex(real64), parameter :: zero = (0.0d0, 0.0d0)

        ! Set f_ to zeros
        f_ = zero

        ! Assign the forcing term to the appropriate node
        f_(5) = applied_force
    end subroutine

! ------------------------------------------------------------------------------
    subroutine plot_deflection(nodes_, u_, leg_nodes_, scaling, title)
        use fplot_core
        !! Plots the deformed shape
        class(node), intent(in), dimension(:) :: nodes_
            !! The node list.
        real(real64), intent(in), dimension(:) :: u_
            !! The nodal deformation vector.
        integer(int32), intent(in), dimension(:,:) :: leg_nodes_
            !! The global node index of each node along each leg, in order.
        real(real64), intent(in) :: scaling
            !! A scaling factor to aid in visualization.
        character(len = *), intent(in) :: title
            !! The plot title.

        ! Local Variables
        integer(int32) :: j, leg, nid, npts
        real(real64), allocatable, dimension(:) :: x, y, dx, dy, umag
        type(plot_2d) :: plt
        type(plot_data_2d), allocatable, dimension(:) :: pd
        type(rainbow_colormap) :: map

        ! Initialization
        npts = size(leg_nodes_, 2)
        allocate(pd(size(leg_nodes_, 1)))
        allocate(x(npts), y(npts), dx(npts), dy(npts), umag(npts))

        call plt%initialize()
        call plt%set_colormap(map)
        call plt%set_title(title)

        ! Plot each leg as a single line, undeformed and deformed
        do leg = 1, size(leg_nodes_, 1)
            do j = 1, npts
                nid = leg_nodes_(leg, j)
                x(j) = nodes_(nid)%x
                y(j) = nodes_(nid)%y
                dx(j) = u_((nid - 1) * 3 + 1)
                dy(j) = u_((nid - 1) * 3 + 2)
            end do
            umag = sqrt(dx**2 + dy**2)

            call plt%push(x, y, lw = 1.0, lc = CLR_BLACK, ls = LINE_DOTTED)

            call pd(leg)%define_data(x + scaling * dx, y + scaling * dy, &
                c = umag)
            call pd(leg)%set_line_width(4.0)
            call plt%push(pd(leg))
        end do

        call plt%draw()
    end subroutine

! ------------------------------------------------------------------------------
    subroutine plot_frf(rsp_, dof_)
        use fplot_core
        !! Plots the FRF of the requested degree of freedom.
        type(frf), intent(in) :: rsp_
            !! The frequency response.
        integer(int32), intent(in) :: dof_
            !! The degree of freedom to plot.

        ! Local Variables
        type(multiplot) :: plt
        type(plot_2d) :: plt1, plt2

        ! Create the plot
        call plt%initialize(2, 1, width = 1000, height = 500)
        call plt1%initialize()
        call plt2%initialize()

        call plt1%set_x_axis_title("f [Hz]")
        call plt1%set_y_axis_title("|X| [dB]")
        call plt2%set_x_axis_title("f [Hz]")
        call plt2%set_y_axis_title("{/Symbol f} [deg]")

        call plt1%push( &
            rsp_%frequency / (2.0d0 * pi), &
            2.0d1 * log10(abs(rsp_%responses(:,dof_)) / abs(rsp_%responses(1,dof_))) &
        )
        call plt%set(1, 1, plt1)

        call plt2%push( &
            rsp_%frequency / (2.0d0 * pi), &
            (1.8d2 / pi) * atan2(aimag(rsp_%responses(:,dof_)), real(rsp_%responses(:,dof_))) &
        )
        call plt%set(2, 1, plt2)
        
        call plt%draw()
    end subroutine

! ------------------------------------------------------------------------------
end program
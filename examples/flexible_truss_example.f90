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
        flexible_material, linear_elastic_isotropic_material, solve_static_system, &
        apply_force_boundary_conditions, restore_constrained_values, degree_of_freedom
    use fplot_core, only : plot_2d, plot_data_2d, plot_3d, plot_data_3d, plot_axis, rainbow_colormap, CLR_BLACK, LINE_DOTTED
    use linalg, only : csr_matrix, size, eigen, matmul
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
    ! DOFs: three translations and three rotations. Both analyses use the
    ! tangent stiffness at the undeformed configuration.

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
    real(real64), parameter :: applied_force = 1.0d4
    integer, parameter :: n_plot_modes = 6
    real(real64), parameter :: mode_plot_amplitude = 0.05d0

    type(linear_elastic_isotropic_material), target :: material
    class(flexible_material), pointer :: mat
    type(finite_element_mesh) :: mesh
    type(csr_matrix) :: mass
    type(csr_matrix) :: stiffness
    integer, allocatable, dimension(:) :: bc
    real(real64), allocatable, dimension(:) :: q
    real(real64), allocatable, dimension(:) :: force
    real(real64), allocatable, dimension(:) :: u
    real(real64), allocatable, dimension(:) :: freqs
    real(real64), allocatable, dimension(:,:) :: modes
    type(degree_of_freedom), allocatable, dimension(:) :: translations
    integer, dimension(3, nnodes) :: translation_eq
    integer, dimension(nel + 1, nlegs) :: leg_nodes
    integer :: ndof
    integer :: i
    integer :: leg
    integer :: el
    integer :: node_count
    integer :: elem_count
    real(real64) :: t
    real(real64) :: amplitude
    real(real64), dimension(3) :: start
    real(real64), dimension(3) :: finish
    character(len = 64) :: title

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

    ! Apply the apex load in global -Y, then remove the same support equations
    ! from the force vector and both matrices.
    call apply_force_boundary_conditions(mesh, [nodal_condition(2, DOF_UY, -applied_force)], force)
    force = apply_boundary_conditions(bc, force)
    mass = apply_boundary_conditions(bc, mass)
    stiffness = apply_boundary_conditions(bc, stiffness)

    ! Resolve global translations once for plotting the full solution vectors.
    call mesh%resolve_dofs([(i, i=1, nnodes), (i, i=1, nnodes), (i, i=1, nnodes)], &
        [(DOF_UX, i=1, nnodes), (DOF_UY, i=1, nnodes), (DOF_UZ, i=1, nnodes)], translations)
    translation_eq(1,:) = translations(:nnodes)%equation_id
    translation_eq(2,:) = translations(nnodes+1:2*nnodes)%equation_id
    translation_eq(3,:) = translations(2*nnodes+1:)%equation_id

    print "(A,I0,A,I0,A,I0,A)", "Mesh: ", size(mesh%nodes), " nodes, ", size(mesh%elements), &
        " elements, ", size(q), " DOFs"
    print "(A,I0)", "Free DOFs after applying supports: ", size(stiffness, 1)
    
    ! Restore zero translations and rotations at the supports before plotting.
    q = solve_static_system(stiffness, force)
    u = restore_constrained_values(bc, q)
    print "(A,ES13.5,A)", "Static apex UY: ", u(translation_eq(2,2)), " m"

    ! Solve K*phi = omega**2*M*phi using the same reduced reference tangent.
    call solve_modes()
    freqs = freqs / (2.0d0 * pi)
    do i = 1, size(freqs)
        print "(A,I0,A,F0.3,A)", "Mode ", i, ": ", freqs(i), " Hz"
    end do

    call plot_mesh()
    call plot_deflection(u, 1.0d2, "Static Deformation (x100)")
    do i = 1, size(freqs)
        u = restore_constrained_values(bc, modes(:,i))
        ! Eigenvector amplitudes are arbitrary. Scale the largest translation
        ! to a fixed display amplitude, without changing the solved mode.
        amplitude = maxval(sqrt(u(translation_eq(1,:))**2 + u(translation_eq(2,:))**2 + u(translation_eq(3,:))**2))
        if (amplitude <= 0.0d0) error stop "Mode has no translational displacement to plot"
        write(title, "(A,I0,A,F0.3,A)") "Mode ", i, ": ", freqs(i), " Hz"
        call plot_mode_shape(u, mode_plot_amplitude / amplitude, trim(title))
    end do

contains
! ------------------------------------------------------------------------------
    subroutine solve_modes()
        !! Extract the six lowest modes of the sparse tangent/mass pencil.
        real(real64), parameter :: sigma = -(2.0d0 * pi * 1000.0d0)**2
        real(real64), parameter :: residual_tolerance = 1.0d-4
        type(csr_matrix) :: scaled_mass
        type(csr_matrix) :: scaled_stiffness
        real(real64), dimension(size(stiffness, 1)) :: scale
        real(real64), allocatable, dimension(:) :: eigenvalues
        real(real64), allocatable, dimension(:,:) :: eigenvectors
        real(real64), dimension(size(stiffness, 1)) :: k_phi
        real(real64), dimension(size(stiffness, 1)) :: lambda_m_phi
        integer :: mode
        integer :: k
        real(real64) :: residual

        ! Equilibrate the stiffness diagonal with a congruence transform.
        ! This preserves the spectrum and improves the shift-invert solves.
        do k = 1, size(scale)
            scale(k) = stiffness%get(k, k)
        end do
        if (any(scale <= 0.0d0)) error stop "Nonpositive tangent stiffness diagonal"
        scale = 1.0d0 / sqrt(scale)
        scaled_mass = scale_matrix(mass, scale)
        scaled_stiffness = scale_matrix(stiffness, scale)

        ! A negative shift improves conditioning without changing which
        ! positive eigenvalues are closest to the shift.
        call eigen(scaled_stiffness, scaled_mass, n_plot_modes, eigenvalues, eigenvectors, sigma=sigma)
        if (size(eigenvalues) /= n_plot_modes .or. any(eigenvalues <= 0.0d0)) &
            error stop "Expected six positive modal eigenvalues"
        allocate(freqs(n_plot_modes), modes(size(stiffness, 1), n_plot_modes))
        do mode = 1, n_plot_modes
            k = minloc(eigenvalues, dim=1)
            freqs(mode) = sqrt(eigenvalues(k))
            modes(:,mode) = scale * eigenvectors(:,k)
            k_phi = matmul(stiffness, modes(:,mode))
            lambda_m_phi = eigenvalues(k) * matmul(mass, modes(:,mode))
            residual = norm2(k_phi - lambda_m_phi) / norm2(lambda_m_phi)
            if (.not.(residual <= residual_tolerance)) then
                print "(A,I0,A,ES13.5)", "Mode ", mode, " relative residual: ", residual
                error stop "Modal eigenpair residual exceeds tolerance"
            end if
            eigenvalues(k) = huge(0.0d0)
        end do
    end subroutine

! ------------------------------------------------------------------------------
    function scale_matrix(matrix, scale) result(scaled)
        type(csr_matrix), intent(in) :: matrix
        real(real64), intent(in), dimension(:) :: scale
        type(csr_matrix) :: scaled
        integer :: row
        integer :: k

        scaled = matrix
        do row = 1, size(scale)
            do k = matrix%row_indices(row), matrix%row_indices(row+1) - 1
                scaled%values(k) = matrix%values(k) * (scale(row) * scale(matrix%column_indices(k)))
            end do
        end do
    end function

! ------------------------------------------------------------------------------
    subroutine plot_deflection(displacement, scaling, title)
        !! Overlay the reference centerlines and displacement-colored shapes.
        real(real64), intent(in), dimension(:) :: displacement
        real(real64), intent(in) :: scaling
        character(len = *), intent(in) :: title
        type(plot_2d) :: plt
        type(plot_data_2d), dimension(nlegs) :: data
        type(rainbow_colormap) :: map
        real(real64), dimension(3, nel + 1) :: xyz
        real(real64), dimension(3, nel + 1) :: d
        integer :: leg

        call plt%initialize()
        call plt%set_title(title)
        call plt%set_colormap(map)
        call plt%set_x_axis_title("X (m)")
        call plt%set_y_axis_title("Y (m)")
        call plt%show_legend(.false.)

        do leg = 1, nlegs
            call leg_displacements(leg, displacement, xyz, d)
            call plt%push(xyz(1,:), xyz(2,:), lw=1.0, lc=CLR_BLACK, ls=LINE_DOTTED)
            call data(leg)%define_data(xyz(1,:) + scaling*d(1,:), xyz(2,:) + scaling*d(2,:), c=norm2(d, dim=1))
            call data(leg)%set_line_width(4.0)
            call plt%push(data(leg))
        end do
        call plt%draw()
    end subroutine

! ------------------------------------------------------------------------------
    subroutine plot_mode_shape(displacement, scaling, title)
        !! Overlay the mode on the undeformed truss in an upright 3D view.
        real(real64), intent(in), dimension(:) :: displacement
        real(real64), intent(in) :: scaling
        character(len = *), intent(in) :: title
        type(plot_3d) :: plt
        type(plot_data_3d), dimension(nlegs) :: data
        class(plot_axis), pointer :: depth_axis
        type(rainbow_colormap) :: map
        real(real64), dimension(3, nel + 1) :: xyz
        real(real64), dimension(3, nel + 1) :: d
        real(real64), dimension(3, nel + 1) :: deformed
        integer :: leg

        call plt%initialize()
        call plt%set_title(title)
        call plt%set_colormap(map)
        call plt%set_x_axis_title("X (m)")
        call plt%set_y_axis_title("-Z (m)")
        call plt%set_z_axis_title("Y (m)")
        call plt%set_elevation(60.0d0)
        call plt%set_azimuth(30.0d0)
        call plt%show_legend(.false.)

        ! Use a common depth range so autoscaling doesn't magnify numerical noise
        ! in the out-of-plane component of an otherwise planar mode.
        depth_axis => plt%get_y_axis()
        call depth_axis%set_limits(-mode_plot_amplitude, mode_plot_amplitude)
        do leg = 1, nlegs
            call leg_displacements(leg, displacement, xyz, d)
            deformed = xyz + scaling*d
            ! Rotate display coordinates +90 degrees about X: (X, Y, Z) -> (X, -Z, Y).
            call plt%push(xyz(1,:), -xyz(3,:), xyz(2,:), lw=1.0, lc=CLR_BLACK, ls=LINE_DOTTED)
            call data(leg)%define_data(deformed(1,:), -deformed(3,:), deformed(2,:), c=norm2(d, dim=1))
            call data(leg)%set_line_width(4.0)
            call plt%push(data(leg))
        end do
        call plt%draw()
    end subroutine

! ------------------------------------------------------------------------------
    subroutine leg_displacements(leg, displacement, xyz, d)
        !! Gather reference coordinates and global translations for one leg.
        integer, intent(in) :: leg
        real(real64), intent(in), dimension(:) :: displacement
        real(real64), intent(out), dimension(3, nel + 1) :: xyz
        real(real64), intent(out), dimension(3, nel + 1) :: d
        integer :: j
        integer :: nid

        do j = 1, nel + 1
            nid = leg_nodes(j,leg)
            xyz(:,j) = mesh%nodes(nid)%coordinates
            d(:,j) = displacement(translation_eq(:,nid))
        end do
    end subroutine

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
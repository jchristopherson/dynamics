program example
    use iso_fortran_env
    use dynamics
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
    ! acting downwards.  Each leg of the structure will be constructed of a
    ! single 2D beam element resulting in a total of 5 elements.  Because of the
    ! small size of this problem, dense matrices will be utilized.
    !
    ! Each node has 3 degrees-of-freedom: the two in-plane translations along 
    ! with one rotation; therefore, with 4 nodes, the whole system has 12 
    ! degrees-of-freedom.

    ! Constants
    character, parameter :: tab = achar(9)
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

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

    ! Variables
    type(material) :: mat
    type(node) :: nodes(4)
    type(beam_element_2d) :: beams(5)
    integer(int32) :: i, bc(4)
    real(real64), allocatable, dimension(:) :: F, Fbc, ubc, u, Fr, freqs
    real(real64), allocatable, dimension(:,:) :: K, Kbc, M, Mbc, shapesbc, shapes
    character(len = 32) :: str

    ! Create the material
    mat = material(modulus, poissons_ratio, density)

    ! Create the list of nodes
    nodes(1) = node(1, 3, 0.0d0, 0.0d0, 0.0d0)
    nodes(2) = node(2, 3, 0.5d0, 0.2d0, 0.0d0)
    nodes(3) = node(3, 3, 1.0d0, 0.0d0, 0.0d0)
    nodes(4) = node(4, 3, 0.5d0, 0.0d0, 0.0d0)

    ! Create the element list
    beams(1) = beam_element_2d(mat, area, moi, nodes(1), nodes(2))
    beams(2) = beam_element_2d(mat, area, moi, nodes(2), nodes(3))
    beams(3) = beam_element_2d(mat, area, moi, nodes(3), nodes(4))
    beams(4) = beam_element_2d(mat, area, moi, nodes(4), nodes(1))
    beams(5) = beam_element_2d(mat, area, moi, nodes(2), nodes(4))

    ! Assemble the stiffness and mass matrices
    call assemble_dynamic_system(12, beams, nodes, M, K)

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
    call plot_deflection(nodes, u, 1.0d2, "Static Deformation")

! ******************************************************************************
! EIGEN ANALYSIS
! ------------------------------------------------------------------------------
    ! Compute the resonant modes of the system
    call modal_response(Mbc, Kbc, freqs, shapesbc)

    ! Convert the frequency units to units of Hz, from rad/s
    freqs = freqs / (2.0d0 * pi)

    ! Reconstruct each mode shape, including the boundary conditions
    allocate(shapes(size(u), size(shapesbc, 2)))
    do i = 1, size(shapesbc, 2)
        shapes(:,i) = restore_constrained_values(bc, shapesbc(:,i))
    end do

    ! Write out each mode shape
    print "(A)", new_line('a') // "MODAL RESPONSE"
    do i = 1, size(freqs)
        print "(A, I0, A, F0.3, A)", "Mode ", i, ": ", freqs(i), " Hz"
    end do

    ! Plot each mode shape
    do i = 1, size(freqs)
        write(str, "(A, I0, A, F0.3, A)") "Mode ", i, ": ", freqs(i), " Hz"
        call plot_deflection(nodes, shapes(:,i), 1.0d0, trim(str))
    end do

contains
    subroutine plot_deflection(nodes_, u_, scaling, title)
        use fplot_core
        !! Plots the deformed shape
        class(node), intent(in), dimension(:) :: nodes_
            !! The node list.
        real(real64), intent(in), dimension(:) :: u_
            !! The nodal deformation vector.
        real(real64), intent(in) :: scaling
            !! A scaling factor to aid in visualization.
        character(len = *), intent(in) :: title
            !! The plot title.

        ! Local Variables
        real(real64) :: &
            x1(2), x2(2), x3(2), x4(2), x5(2), &
            y1(2), y2(2), y3(2), y4(2), y5(2), &
            dx1(2), dx2(2), dx3(2), dx4(2), dx5(2), &
            dy1(2), dy2(2), dy3(2), dy4(2), dy5(2), &
            u1(2), u2(2), u3(2), u4(2), u5(2)
        type(plot_2d) :: plt
        type(plot_data_2d) :: pd1, pd2, pd3, pd4, pd5
        type(rainbow_colormap) :: map

        ! Generate the deformation vectors
        dx1 = [u_(1), u_(4)]
        dy1 = [u_(2), u_(5)]
        dx2 = [u_(4), u_(7)]
        dy2 = [u_(5), u_(8)]
        dx3 = [u_(7), u_(10)]
        dy3 = [u_(8), u_(11)]
        dx4 = [u_(10), u_(1)]
        dy4 = [u_(11), u_(2)]
        dx5 = [u_(10), u_(4)]
        dy5 = [u_(11), u_(5)]

        u1 = sqrt(dx1**2 + dy1**2)
        u2 = sqrt(dx2**2 + dy2**2)
        u3 = sqrt(dx3**2 + dy3**2)
        u4 = sqrt(dx4**2 + dy4**2)
        u5 = sqrt(dx5**2 + dy5**2)

        ! Generate the mesh lines
        x1 = [nodes_(1)%x, nodes_(2)%x]
        y1 = [nodes_(1)%y, nodes_(2)%y]
        x2 = [nodes_(2)%x, nodes_(3)%x]
        y2 = [nodes_(2)%y, nodes_(3)%y]
        x3 = [nodes_(3)%x, nodes_(4)%x]
        y3 = [nodes_(3)%y, nodes_(4)%y]
        x4 = [nodes_(4)%x, nodes_(1)%x]
        y4 = [nodes_(4)%y, nodes_(1)%y]
        x5 = [nodes_(4)%x, nodes_(2)%x]
        y5 = [nodes_(4)%y, nodes_(2)%y]

        ! Generate the mesh plot
        call plt%initialize()
        call plt%set_colormap(map)
        call plt%set_title(title)
        call plt%push(x1, y1, lw = 1.0, lc = CLR_BLACK, ls = LINE_DOTTED)
        call plt%push(x2, y2, lw = 1.0, lc = CLR_BLACK, ls = LINE_DOTTED)
        call plt%push(x3, y3, lw = 1.0, lc = CLR_BLACK, ls = LINE_DOTTED)
        call plt%push(x4, y4, lw = 1.0, lc = CLR_BLACK, ls = LINE_DOTTED)
        call plt%push(x5, y5, lw = 1.0, lc = CLR_BLACK, ls = LINE_DOTTED)

        call pd1%define_data(x1 + scaling * dx1, y1 + scaling * dy1, c = u1)
        call pd1%set_line_width(4.0)
        call plt%push(pd1)

        call pd2%define_data(x2 + scaling * dx2, y2 + scaling * dy2, c = u2)
        call pd2%set_line_width(4.0)
        call plt%push(pd2)

        call pd3%define_data(x3 + scaling * dx3, y3 + scaling * dy3, c = u3)
        call pd3%set_line_width(4.0)
        call plt%push(pd3)

        call pd4%define_data(x4 + scaling * dx4, y4 + scaling * dy4, c = u4)
        call pd4%set_line_width(4.0)
        call plt%push(pd4)

        call pd5%define_data(x5 + scaling * dx5, y5 + scaling * dy5, c = u5)
        call pd5%set_line_width(4.0)
        call plt%push(pd5)

        call plt%draw()
    end subroutine
end program
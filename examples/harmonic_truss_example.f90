program harmonic_truss_example
    use iso_fortran_env, only : int32, real64
    use dynamics, only : dense_generalized_alpha_integrator, apply_boundary_conditions
    use fplot_core
    implicit none

    ! A pin-jointed, five-bar planar truss (all lengths in metres):
    !       2
    !      /|\
    !     / | \
    !    1--4--3
    ! Node 1 is pinned (x, y); node 3 is supported vertically. A vertical
    ! harmonic force acts at the apex, node 2. Bars carry axial force only.
    integer(int32), parameter :: nnode = 4, ndof = 2 * nnode, nbar = 5
    integer(int32), parameter :: nsteps = 2000
    real(real64), parameter :: pi = acos(-1.0d0)
    real(real64), parameter :: dt = 2.0d-4, frequency = 20.0d0
    real(real64), parameter :: force_amplitude = 2.0d3
    real(real64), parameter :: elastic_modulus = 70.0d9, density = 2.7d3
    real(real64), parameter :: area = 2.5d-4
    real(real64), parameter :: positions(2,nnode) = reshape([ &
        0.0d0, 0.0d0,  0.5d0, 0.3d0,  1.0d0, 0.0d0,  0.5d0, 0.0d0], [2,nnode])
    integer(int32), parameter :: ends(2,nbar) = reshape([ &
        1, 2,  2, 3,  3, 4,  4, 1,  2, 4], [2,nbar])

    integer(int32) :: bar, time_index, apex_y, middle_y, middle_x, right_x
    integer(int32) :: restrained(3)
    real(real64) :: mass(ndof,ndof), stiffness(ndof,ndof)
    real(real64), allocatable, dimension(:,:) :: reduced_mass, reduced_stiffness, damping
    real(real64), allocatable, dimension(:) :: displacement, velocity, acceleration
    real(real64), allocatable, dimension(:) :: force_current, force_next
    real(real64) :: time(nsteps+1), apex_motion(nsteps+1), middle_motion(nsteps+1)
    real(real64) :: middle_horizontal(nsteps+1), right_horizontal(nsteps+1)
    type(dense_generalized_alpha_integrator) :: integrator
    type(multiplot) :: plt
    type(plot_2d) :: plt1, plt2
    class(legend), pointer :: lgnd

    mass = 0.0d0
    stiffness = 0.0d0
    do bar = 1, nbar
        call assemble_bar(ends(1,bar), ends(2,bar), positions, &
            elastic_modulus, density, area, mass, stiffness)
    end do

    restrained = [1, 2, 6]
    reduced_mass = apply_boundary_conditions(restrained, mass)
    reduced_stiffness = apply_boundary_conditions(restrained, stiffness)
    damping = 1.0d-2 * reduced_mass + 2.0d-5 * reduced_stiffness

    ! Each global node has x and y DOFs. Removing supports shifts their
    ! indices in the reduced vectors but preserves the order of free DOFs.
    apex_y = 4 - count(restrained < 4)
    middle_x = 7 - count(restrained < 7)
    middle_y = 8 - count(restrained < 8)
    right_x = 5 - count(restrained < 5)
    allocate(displacement(size(reduced_mass,1)), velocity(size(reduced_mass,1)), &
        acceleration(size(reduced_mass,1)), force_current(size(reduced_mass,1)), &
        force_next(size(reduced_mass,1)), source = 0.0d0)

    ! The initial load and state are zero, so the initial acceleration is
    ! consistent with M*a + C*v + K*u = f at t = 0.
    call integrator%initialize(reduced_mass, damping, reduced_stiffness, rho_infinity = 0.7d0)
    do time_index = 1, nsteps + 1
        time(time_index) = (time_index - 1) * dt
        apex_motion(time_index) = 1.0d3 * displacement(apex_y)
        middle_motion(time_index) = 1.0d3 * displacement(middle_y)
        middle_horizontal(time_index) = 1.0d3 * displacement(middle_x)
        right_horizontal(time_index) = 1.0d3 * displacement(right_x)
        if (time_index > nsteps) exit

        force_next = 0.0d0
        force_next(apex_y) = -force_amplitude * &
            sin(2.0d0 * pi * frequency * real(time_index, real64) * dt)
        call integrator%step(force_current, force_next, dt, &
            displacement, velocity, acceleration)
        force_current = force_next
    end do

    print "(A,F7.2,A)", "Forcing frequency: ", frequency, " Hz"
    print "(A,ES12.4,A)", "Peak apex vertical motion: ", maxval(abs(apex_motion)), " mm"
    print "(A,ES12.4,A)", "Peak midspan vertical motion: ", maxval(abs(middle_motion)), " mm"

    call plt%initialize(2, 1, width = 1100, height = 750)
    call plt1%initialize()
    call plt2%initialize()
    lgnd => plt2%get_legend()
    call lgnd%set_is_visible(.true.)

    call plt1%set_title("Vertical displacement")
    call plt1%set_x_axis_title("Time [s]")
    call plt1%set_y_axis_title("Displacement [mm]")
    call plt1%push(time, apex_motion, name = "Apex (node 2)", lw = 2.0)
    call plt1%push(time, middle_motion, name = "Midspan (node 4)", lw = 2.0)
    call plt%set(1, 1, plt1)

    call plt2%set_title("Horizontal displacement")
    call plt2%set_x_axis_title("Time [s]")
    call plt2%set_y_axis_title("Displacement [mm]")
    call plt2%push(time, middle_horizontal, name = "Midspan (node 4)", lw = 2.0)
    call plt2%push(time, right_horizontal, name = "Roller (node 3)", lw = 2.0)
    call plt%set(2, 1, plt2)

    call plt%draw()

contains

subroutine assemble_bar(first, second, xy, modulus, material_density, cross_section, m, k)
    integer(int32), intent(in) :: first, second
    real(real64), intent(in), dimension(:,:) :: xy
    real(real64), intent(in) :: modulus, material_density, cross_section
    real(real64), intent(inout), dimension(:,:) :: m, k
    integer(int32) :: local_row, local_col, global_dof(4)
    real(real64) :: length, direction(2), axial(4), element_mass

    direction = xy(:,second) - xy(:,first)
    length = norm2(direction)
    direction = direction / length
    axial = [-direction(1), -direction(2), direction(1), direction(2)]
    global_dof = [2*first-1, 2*first, 2*second-1, 2*second]
    element_mass = material_density * cross_section * length / 6.0d0

    do local_col = 1, 4
        do local_row = 1, 4
            k(global_dof(local_row),global_dof(local_col)) = &
                k(global_dof(local_row),global_dof(local_col)) + &
                modulus * cross_section / length * axial(local_row) * axial(local_col)
        end do
    end do
    do local_row = 1, 2
        m(global_dof(local_row),global_dof(local_row)) = &
            m(global_dof(local_row),global_dof(local_row)) + 2.0d0 * element_mass
        m(global_dof(local_row+2),global_dof(local_row+2)) = &
            m(global_dof(local_row+2),global_dof(local_row+2)) + 2.0d0 * element_mass
        m(global_dof(local_row),global_dof(local_row+2)) = &
            m(global_dof(local_row),global_dof(local_row+2)) + element_mass
        m(global_dof(local_row+2),global_dof(local_row)) = &
            m(global_dof(local_row+2),global_dof(local_row)) + element_mass
    end do
end subroutine

end program
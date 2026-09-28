program example
    ! This example compares the results for a system of two mass-spring-damper
    ! pairs in series with excitation provided via a force acting directly on
    ! mass 2.  The structural FE approach will be compared directly with a
    ! traditional ODE solver (via the DIFFEQ library).
    use iso_fortran_env
    use diffeq
    use fplot_core
    use dynamics

    ! Constants
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Model Parameters
    real(real64), parameter :: mass_1 = 1.5d0
    real(real64), parameter :: mass_2 = 2.5d0
    real(real64), parameter :: stiffness = 1.5d4
    real(real64), parameter :: damping = 5.0d2
    real(real64), parameter :: frequency = 2.0d1
    real(real64), parameter :: force_amplitude = 2.0d3

    ! Analysis Parameters
    real(real64), parameter :: tmax = 5.0d-1
    real(real64), parameter :: dt = 1.0d-3

    ! Structural Variables
    type(node) :: nodes(3)
    type(spring_element_2d) :: springs(2)
    type(damper_element_2d) :: dampers(2)
    type(mass_element_2d) :: masses(2)
    integer(int32) :: i, n, bc(4)
    real(real64) :: M(2, 2), B(2, 2), K(2, 2), x(2), v(2), a(2), f(2), fnew(2)
    real(real64), allocatable, dimension(:) :: t
    real(real64), allocatable, dimension(:,:) :: M_full, B_full, K_full, &
        position, velocity, acceleration
    type(dense_generalized_alpha_integrator) :: struct_integrator

    ! ODE Solver Variables
    type(ode_container) :: model
    type(runge_kutta_45) :: integrator
    real(real64), allocatable, dimension(:,:) :: solution

    ! Plot Variables
    type(multiplot) :: plt
    type(plot_2d) :: plt1, plt2, plt3, plt4
    class(legend), pointer :: lgnd

    ! Define the nodal locations
    nodes(1) = node(1, 2, 0.0d0, 0.0d0, 0.0d0)
    nodes(2) = node(2, 2, 0.0d0, 1.0d0, 0.0d0)
    nodes(3) = node(3, 2, 0.0d0, 2.0d0, 0.0d0)

    ! Define the spring elements
    springs(1) = spring_element_2d( &
        stiffness, &
        nodes(1), &
        nodes(2) &
    )
    springs(2) = spring_element_2d( &
        stiffness, &
        nodes(2), &
        nodes(3) &
    )

    ! Define the damper elements
    dampers(1) = damper_element_2d( &
        damping, &
        nodes(1), &
        nodes(2) &
    )
    dampers(2) = damper_element_2d( &
        damping, &
        nodes(2), &
        nodes(3) &
    )

    ! Define the mass elements
    masses(1) = mass_element_2d( &
        mass_1, &
        nodes(2) &
    )
    masses(2) = mass_element_2d( &
        mass_2, &
        nodes(3) &
    )

    ! Assemble the system matrices
    call assemble_discrete_system(6, masses, dampers, springs, nodes, &
        M_full, B_full, K_full)

    ! Apply boundary conditions (fix node 1 and only allow y-direction 
    ! translations on the other nodes).
    bc = [1, 2, 3, 5]
    M = apply_boundary_conditions(bc, M_full)
    B = apply_boundary_conditions(bc, B_full)
    K = apply_boundary_conditions(bc, K_full)

    ! Define a time vector
    n = floor(tmax / dt) + 1
    t = [(i * dt, i = 0, n - 1)]

    ! Set up the integrator
    call struct_integrator%initialize(M, B, K, rho_infinity = 7.0d-1)

    ! Compute the solution for the structural model - assume zero initial
    ! conditions
    allocate(position(n, 2), velocity(n, 2), acceleration(n, 2), source = 0.0d0)
    x = 0.0d0
    v = 0.0d0
    a = 0.0d0
    call force_vector(0.0d0, f)
    do i = 2, n
        ! Compute the force vector at t(n+1)
        call force_vector(t(i), fnew)

        ! Compute the solution at t(n+1)
        call struct_integrator%step(f, fnew, dt, x, v, a)

        ! Store the results
        position(i,:) = x
        velocity(i,:) = v
        acceleration(i,:) = a

        ! Update the forcing term
        f = fnew
    end do

! ******************************************************************************
! ODE Integrator
    model%fcn => equations
    call integrator%solve(model, [0.0d0, tmax], [0.0d0, 0.0d0, 0.0d0, 0.0d0])
    solution = integrator%get_solution()

! ******************************************************************************
! Plotting
    call plt%initialize(2, 2, width = 1200, height = 800)
    call plt1%initialize()
    call plt2%initialize()
    call plt3%initialize()
    call plt4%initialize()
    lgnd => plt1%get_legend()
    call lgnd%set_is_visible(.true.)
    call lgnd%set_draw_border(.false.)

    call plt1%push(t, 1.0d3 * position(:,1), name = "Structural")
    call plt1%push(solution(:,1), 1.0d3 * solution(:,2), name = "ODE", ls = LINE_DASHED)
    call plt1%set_x_axis_title("t [sec]")
    call plt1%set_y_axis_title("x_1 [mm]")
    call plt%set(1, 1, plt1)

    call plt2%push(t, 1.0d3 * position(:,2), name = "Structural")
    call plt2%push(solution(:,1), 1.0d3 * solution(:,4), name = "ODE", ls = LINE_DASHED)
    call plt2%set_x_axis_title("t [sec]")
    call plt2%set_y_axis_title("x_2 [mm]")
    call plt%set(2, 1, plt2)

    call plt3%push(t, velocity(:,1), name = "Structural")
    call plt3%push(solution(:,1), solution(:,3), name = "ODE", ls = LINE_DASHED)
    call plt3%set_x_axis_title("t [sec]")
    call plt3%set_y_axis_title("v_1 [m/s]")
    call plt%set(1, 2, plt3)

    call plt4%push(t, velocity(:,2), name = "Structural")
    call plt4%push(solution(:,1), solution(:,5), name = "ODE", ls = LINE_DASHED)
    call plt4%set_x_axis_title("t [sec]")
    call plt4%set_y_axis_title("v_2 [m/s]")
    call plt%set(2, 2, plt4)

    call plt%draw()
contains
    subroutine force_vector(t_, f_)
        !! A routine for computing the external forcing function.
        real(real64), intent(in) :: t_
            !! The current simulation time.
        real(real64), intent(out), dimension(:) :: f_
            !! The external force vector.

        f_(1) = 0.0d0
        f_(2) = force_amplitude * sin(2.0d0 * pi * frequency * t_)
    end subroutine

    subroutine equations(t_, state_, derivatives_, args_)
        !! The equations to integrate via a traditional ODE integrator.
        real(real64), intent(in) :: t_
            !! The current simulation time.
        real(real64), intent(in), dimension(:) :: state_
            !! The state vector.
        real(real64), intent(out), dimension(:) :: derivatives_
            !! The derivative vector.
        class(*), intent(inout), optional :: args_
            !! User communication container.

        ! Local Variables
        real(real64) :: fvec_(2)

        ! Compute the external force
        call force_vector(t_, fvec_)

        ! Define the derivatives
        derivatives_(1) = state_(2)
        derivatives_(2) = (fvec_(1) - ( &
            2.0d0 * damping * state_(2) - damping * state_(4) + &
            2.0d0 * stiffness * state_(1) - stiffness * state_(3))) / mass_1
        derivatives_(3) = state_(4)
        derivatives_(4) = (fvec_(2) - ( &
            -damping * state_(2) + damping * state_(4) - &
            stiffness * state_(1) + stiffness * state_(3))) / mass_2
    end subroutine
end program
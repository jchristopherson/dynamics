module equation_container
    use iso_fortran_env
    implicit none

    ! Parameters
    real(real64), parameter :: sigma = 10.0d0
    real(real64), parameter :: rho = 28.0d0
    real(real64), parameter :: beta = 8.0d0 / 3.0d0

contains
    ! Lorenz Equation
    subroutine eom(t, x, dxdt, args)
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: dxdt
        class(*), intent(inout), optional :: args

        ! Model
        dxdt(1) = sigma * (x(2) - x(1))
        dxdt(2) = x(1) * (rho - x(3)) - x(2)
        dxdt(3) = x(1) * x(2) - beta * x(3)
    end subroutine
end module

program main
    use iso_fortran_env
    use equation_container
    use diffeq
    use dynamics_maps
    use dynamics_geometry
    use fplot_core
    implicit none

    ! Local Variables
    type(ode_container) :: mdl
    type(runge_kutta_45) :: integrator
    real(real64), allocatable, dimension(:,:) :: sol, map
    type(plane) :: pln

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_3d) :: plt3
    type(plot_data_2d) :: pd
    type(plot_data_3d) :: pd3

    ! Define the slice plane (x-y plane offset 25 units in z)
    pln = plane([0.0d0, 0.0d0, 25.0d0], [0.0d0, 0.0d0, 1.0d0]) ! point and normal

    ! Solve the ODE system
    mdl%fcn => eom
    call integrator%solve(mdl, [0.0d0, 1.0d2], [1.0d0, 0.0d0, 1.0d0])
    sol = integrator%get_solution()

    ! Construct the section
    map = poincare_map(sol(:,2), sol(:,3), sol(:,4), pln)

    ! Plot
    call plt%initialize()
    call plt3%initialize()

    ! Trajectory Plot
    call pd3%define_data(sol(:,2), sol(:,3), sol(:,4))
    call plt3%push(pd3)
    call plt3%draw()

    ! Poincare Plot
    call pd%define_data(map(:,1), map(:,2))
    call pd%set_draw_line(.false.)
    call pd%set_draw_markers(.true.)
    call plt%push(pd)
    call plt%draw()

    ! Plot the Poincare map on the original trajectory
    call pd3%define_data(map(:,1), map(:,2), map(:,3))
    call pd3%set_draw_line(.false.)
    call pd3%set_draw_markers(.true.)
    call pd3%set_line_color(CLR_RED)
    call plt3%push(pd3)
    call plt3%draw()
end program
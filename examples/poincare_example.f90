module equation_container
    use iso_fortran_env
    implicit none

    ! Parameters
    real(real64), parameter :: alpha = 1.0d0
    real(real64), parameter :: beta = 5.0d0
    real(real64), parameter :: delta = 2.0d-2
    real(real64), parameter :: gamma = 8.0d0
    real(real64), parameter :: omega = 0.5d0

contains
    ! Duffings Equation
    pure elemental function duffing_acceleration(t, x, dxdt) result(rst)
        real(real64), intent(in) :: t, x, dxdt
        real(real64) :: rst
        rst = gamma * cos(omega * t) - delta * dxdt - alpha * x - beta * x**3
    end function

    subroutine eom(t, x, dxdt, args)
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: x
        real(real64), intent(out), dimension(:) :: dxdt
        class(*), intent(inout), optional :: args

        ! Model
        dxdt(1) = x(2)
        dxdt(2) = duffing_acceleration(t, x(1), x(2))
    end subroutine
end module

program main
    use iso_fortran_env
    use equation_container
    use diffeq
    use dynamics_maps
    use fplot_core
    implicit none

    ! Local Variables
    type(ode_container) :: mdl
    type(runge_kutta_45) :: integrator
    real(real64), allocatable, dimension(:) :: accel
    real(real64), allocatable, dimension(:,:) :: sol, map

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd

    ! Solve the ODE system
    mdl%fcn => eom
    call integrator%set_step_limit(1000000000)
    call integrator%solve(mdl, [0.0d0, 1.0d5], [0.0d0, 0.0d0])
    sol = integrator%get_solution()
    accel = duffing_acceleration(sol(:,1), sol(:,2), sol(:,3))

    ! Construct the section
    map = poincare_map(sol(:,2), sol(:,3), accel)

    ! Plot
    call plt%initialize()
    call pd%define_data(map(:,1), map(:,2))
    call pd%set_draw_line(.false.)
    call pd%set_draw_markers(.true.)
    call plt%push(pd)
    call plt%draw()
end program
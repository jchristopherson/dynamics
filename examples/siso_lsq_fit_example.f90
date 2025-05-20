module equation_container
    use iso_fortran_env
    use dynamics
    implicit none

contains
    subroutine eom(t, x, dxdt, args)
        real(real64), intent(in) :: t               ! the current time value
        real(real64), intent(in) :: x(:)            ! the current state vector
        real(real64), intent(out) :: dxdt(:)        ! the derivatives
        class(*), intent(inout), optional :: args   ! model information

        ! Local Variables
        real(real64) :: zeta, wn, F

        ! Extract the model information
        select type (args)
        class is (model_information)
            wn = args%model(1)
            zeta = args%model(2)
            F = args%excitation%interpolate_value(t)
        end select

        ! The ODE's:
        ! x" + 2 zeta wn x' + wn**2 x = F(t)
        dxdt(1) = x(2)
        dxdt(2) = F - wn * (2.0d0 * zeta * x(2) + wn * x(1))
    end subroutine
end module

program example
    use iso_fortran_env
    use equation_container
    use dynamics
    use diffeq
    use fstats
    use fplot_core
    implicit none

    ! Parameters
    real(real64), parameter :: fs = 1.0d3
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: zeta = 5.0d-2
    real(real64), parameter :: wn = 3.0d2
    real(real64), parameter :: sigma_pct = 1.0d-1
    real(real64), parameter :: amplitude = 1.0d4

    ! Local Variables
    integer(int32) :: i
    real(real64) :: dt, tmax, p(2), ic(2)
    type(dynamic_system_measurement) :: measurements(1)
    procedure(ode), pointer :: fcn
    type(ode_container) :: mdl
    type(runge_kutta_45) :: integrator
    type(model_information) :: info
    type(linear_interpolator), target :: interp
    real(real64), allocatable, dimension(:,:) :: sol
    type(iteration_controls) :: controls

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2
    class(legend), pointer :: lgnd
    
    ! Generate an initial guess
    p = [2.5d2, 1.0d-1]

    ! Allocate memory for the measurement data we're trying to fit
    allocate( &
        measurements(1)%t(npts), &
        measurements(1)%output(npts), &
        measurements(1)%input(npts) &
    )

    ! Generate a time vector at which to sample the system.
    dt = 1.0d0 / fs
    tmax = dt * (npts - 1.0d0)
    measurements(1)%t = (/ (dt * i, i = 0, npts - 1) /)

    ! Define the forcing function at each time point
    measurements(1)%input = amplitude

    ! Generate the solution for the system
    ic = 0.0d0  ! zero-valued initial conditions
    call interp%initialize(measurements(1)%t, measurements(1)%input)
    info%model = [wn, zeta]
    info%excitation => interp
    mdl%fcn => eom
    call integrator%solve(mdl, measurements(1)%t, ic, args = info)
    sol = integrator%get_solution()
    measurements(1)%output = sol(:,2) + &
        box_muller_sample(0.0d0, 5.0d-3, npts) ! additional noise

    ! This is optional, but is illustrated here to show how to adjust solver
    ! tolerances
    call controls%set_to_default()
    controls%change_in_solution_tolerance = 1.0d-12
    controls%residual_tolerance = 1.0d-8

    ! Set up the problem and solve
    fcn => eom
    call siso_model_fit_least_squares(fcn, measurements, ic, p, &
        controls = controls)

    ! Compare the solution and the actual values
    print "(A)", "NATURAL FREQUENCY TERM:"
    print "(AAF8.3A)", achar(9), "Actual: ", wn, " rad/s"
    print "(AAF8.3A)", achar(9), "Computed: ", p(1), " rad/s"
    print "(AAF8.3A)", achar(9), "Difference: ", p(1) - wn, " rad/s"
    
    print "(A)", "DAMPING TERM:"
    print "(AAF6.3)", achar(9), "Actual: ", zeta
    print "(AAF6.3)", achar(9), "Computed: ", p(2)
    print "(AAF6.3)", achar(9), "Difference: ", p(2) - zeta

    ! Re-evaluate the model with the computed parameters
    info%model = p
    call integrator%clear_buffer()
    call integrator%solve(mdl, measurements(1)%t, ic, args = info)
    sol = integrator%get_solution()
    
    ! Plot the results
    call plt%initialize()
    lgnd => plt%get_legend()
    call lgnd%set_is_visible(.true.)
    call pd1%define_data(measurements(1)%t, measurements(1)%output)
    call pd1%set_name("Data")
    call plt%push(pd1)
    call pd2%define_data(measurements(1)%t, sol(:,2))
    call pd2%set_name("Fit")
    call plt%push(pd2)
    call plt%draw()
end program
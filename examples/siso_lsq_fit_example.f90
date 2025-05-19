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
    use fstats, only : box_muller_sample
    implicit none

    ! Parameters
    real(real64), parameter :: fs = 1.0d3
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: zeta = 2.0d-1
    real(real64), parameter :: wn = 3.0d2
    real(real64), parameter :: sigma_pct = 1.0d-1
    real(real64), parameter :: step_amplitude = 1.0d0

    ! Local Variables
    integer(int32) :: i
    real(real64) :: dt, wng(1), zg(1), p(2), ic(2)
    type(dynamic_system_measurement) :: measurements(1)
    procedure(ode), pointer :: fcn
    
    ! Generate an initial guess randomly from a Gaussian distribution
    wng = box_muller_sample(wn, sigma_pct * wn, 1)
    zg = box_muller_sample(zeta, sigma_pct * zeta, 1)
    p = [wng, zg]

    ! Allocate memory for the measurement data we're trying to fit
    allocate( &
        measurements(1)%t(npts), &
        measurements(1)%output(npts), &
        measurements(1)%input(npts) &
    )

    ! Generate a time vector at which to sample the system.
    dt = 1.0d0 / fs
    measurements(1)%t = (/ (dt * i, i = 0, npts - 1) /)

    ! Define the forcing function at each time point
    measurements(1)%input = step_amplitude

    ! Assume the system is exposed to a step function.  The solution
    ! is then given as follows.
    measurements(1)%output = evaluate_step_response(wn, zeta, step_amplitude, &
        measurements(1)%t)

    ! Set up the problem and solve
    fcn => eom
    ic = 0.0d0  ! zero-valued initial conditions
    call siso_model_fit_least_squares(fcn, measurements, ic, p)

    ! Compare the solution and the actual values
    print "(A)", "NATURAL FREQUENCY TERM:"
    print "(AAF7.3A)", achar(9), "Actual: ", wn, " rad/s"
    print "(AAF7.3A)", achar(9), "Initial Guess: ", wng, " rad/s"
    print "(AAF7.3A)", achar(9), "Computed: ", p(1), " rad/s"
    print "(AAF7.3A)", achar(9), "Difference: ", p(1) - wn, " rad/s"
    
    print "(A)", "DAMPING TERM:"
    print "(AAF5.3)", achar(9), "Actual: ", zeta
    print "(AAF5.3)", achar(9), "Initial Guess: ", zg
    print "(AAF5.3)", achar(9), "Computed: ", p(2)
    print "(AAF6.3)", achar(9), "Difference: ", p(2) - zeta
end program
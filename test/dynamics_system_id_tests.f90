module dynamics_system_id_tests
    use iso_fortran_env
    use fortran_test_helper
    use dynamics
    use diffeq
    use fstats
    implicit none

contains
! ------------------------------------------------------------------------------
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

        ! The ODE:
        ! x" + 2 zeta wn x' + wn**2 x = F(t)
        dxdt(1) = x(2)
        dxdt(2) = F - wn * (2.0d0 * zeta * x(2) + wn * x(1))
    end subroutine

! ------------------------------------------------------------------------------
    subroutine constraints(xg, fg, xc, p, fc, args)
        real(real64), intent(in), dimension(:) :: xg, fg, xc, p
        real(real64), intent(out), dimension(:) :: fc
        class(*), intent(inout), optional :: args

        ! Constraint Equations
        fc(1) = p(1)
        fc(2) = p(2)
    end subroutine

! ------------------------------------------------------------------------------
    function test_siso_model_fit_least_squares() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: tol = 1.0d-1
        real(real64), parameter :: fs = 1.0d3
        integer(int32), parameter :: npts = 1000
        real(real64), parameter :: zeta = 5.0d-2
        real(real64), parameter :: wn = 3.0d2
        real(real64), parameter :: sigma_pct = 1.0d-1
        real(real64), parameter :: amplitude = 1.0d4

        ! Local Variables
        integer(int32) :: i
        real(real64) :: dt, tmax, p(2), ic(2)
        type(dynamic_system_measurement) :: measurements(2)
        procedure(ode), pointer :: fcn
        type(ode_container) :: mdl
        type(runge_kutta_45) :: integrator
        type(model_information) :: info
        type(linear_interpolator), target :: interp, interp2
        real(real64), allocatable, dimension(:,:) :: sol
        type(iteration_controls) :: controls

        ! Initialization
        rst = .true.
        
        ! Generate an initial guess
        p = [2.5d2, 1.0d-1]

! ------------------------------------------------------------------------------
! EXCITATION 1
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
        measurements(1)%output = sol(:,2)

    ! --------------------------------------------------------------------------
    ! EXCITATION 2
        allocate( &
            measurements(2)%t(npts), &
            measurements(2)%output(npts), &
            measurements(2)%input(npts) &
        )

        ! Time Vector
        measurements(2)%t = measurements(1)%t

        ! Forcing Term
        measurements(2)%input = amplitude

        ! Generate the solution for the system
        call interp2%initialize(measurements(2)%t, measurements(2)%input)
        info%excitation => interp2
        call integrator%clear_buffer()
        call integrator%solve(mdl, measurements(2)%t, ic, args = info)
        sol = integrator%get_solution()
        measurements(2)%output = sol(:,2)

    ! --------------------------------------------------------------------------
        ! Set tolerances
        call controls%set_to_default()
        controls%change_in_solution_tolerance = 1.0d-12
        controls%residual_tolerance = 1.0d-8

        ! Set up the problem and solve
        fcn => eom
        call siso_model_fit_least_squares(fcn, measurements, ic, p, &
            controls = controls)

        ! Test the parameters
        if (.not.assert(wn, p(1), tol * wn)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_siso_model_fit_least_squares -1"
        end if

        if (.not.assert(zeta, p(2), tol * zeta)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_siso_model_fit_least_squares -2"
        end if
    end function

! ------------------------------------------------------------------------------
    function test_siso_model_fit_least_squares_constrained() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        real(real64), parameter :: tol = 1.0d-1
        real(real64), parameter :: fs = 1.0d3
        integer(int32), parameter :: npts = 1000
        real(real64), parameter :: zeta = 5.0d-2
        real(real64), parameter :: wn = 3.0d2
        real(real64), parameter :: sigma_pct = 1.0d-1
        real(real64), parameter :: amplitude = 1.0d4

        ! Local Variables
        integer(int32) :: i
        real(real64) :: dt, tmax, p(2), ic(2), tc(2), fc(2)
        type(dynamic_system_measurement) :: measurements(2)
        procedure(ode), pointer :: fcn
        type(ode_container) :: mdl
        type(runge_kutta_45) :: integrator
        type(model_information) :: info
        type(linear_interpolator), target :: interp, interp2
        real(real64), allocatable, dimension(:,:) :: sol
        type(iteration_controls) :: controls
        procedure(constraint_equations), pointer :: cfcn

        ! Initialization
        rst = .true.
        cfcn => constraints
        
        ! Generate an initial guess
        p = [2.5d2, 1.0d-1]

! ------------------------------------------------------------------------------
! EXCITATION 1
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
        measurements(1)%output = sol(:,2)

    ! --------------------------------------------------------------------------
    ! EXCITATION 2
        allocate( &
            measurements(2)%t(npts), &
            measurements(2)%output(npts), &
            measurements(2)%input(npts) &
        )

        ! Time Vector
        measurements(2)%t = measurements(1)%t

        ! Forcing Term
        measurements(2)%input = amplitude

        ! Generate the solution for the system
        call interp2%initialize(measurements(2)%t, measurements(2)%input)
        info%excitation => interp2
        call integrator%clear_buffer()
        call integrator%solve(mdl, measurements(2)%t, ic, args = info)
        sol = integrator%get_solution()
        measurements(2)%output = sol(:,2)

    ! --------------------------------------------------------------------------
        ! Set tolerances
        call controls%set_to_default()
        controls%change_in_solution_tolerance = 1.0d-12
        controls%residual_tolerance = 1.0d-8

        ! Define the constraints
        tc = 0.0d0
        fc = [wn, zeta]

        ! Set up the problem and solve
        fcn => eom
        call siso_model_fit_least_squares(fcn, measurements, ic, p, &
            controls = controls, xc = tc, yc = fc, constraints = cfcn)

        ! Test the parameters
        if (.not.assert(wn, p(1), tol * wn)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_siso_model_fit_least_squares_constrained -1"
        end if

        if (.not.assert(zeta, p(2), tol * zeta)) then
            rst = .false.
            print "(A)", "TEST FAILED: test_siso_model_fit_least_squares_constrained -2"
        end if
    end function

! ------------------------------------------------------------------------------
end module
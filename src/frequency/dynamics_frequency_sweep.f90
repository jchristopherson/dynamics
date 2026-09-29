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
module dynamics_frequency_sweep
    use iso_fortran_env, only : int32, real64
    use diffeq, only : ode_container, ode_integrator
    use dynamics_error_handling
    use dynamics_frequency_response, only : frf
    implicit none
    private
    public :: chirp
    public :: ode_excite, harmonic_ode, ode_integrator, frequency_sweep

    interface
        function ode_excite(t) result(rst)
            !! Defines the interface for a ODE excitation function.
            use iso_fortran_env, only : real64
            real(real64), intent(in) :: t
                !! The value of the independent variable at which to evaluate
                !! the excitation function.
            real(real64) :: rst
                !! The result.
        end function

        subroutine harmonic_ode(freq, t, x, dxdt, args)
            !! Defines a system of ODE's exposed to harmonic excitation.
            use iso_fortran_env, only : real64
            real(real64), intent(in) :: freq
                !! The excitation frequency.
            real(real64), intent(in) :: t
                !! The current time step value.
            real(real64), intent(in), dimension(:) :: x
                !! The value of the solution estimate at time t.
            real(real64), intent(out), dimension(:) :: dxdt
                !! The derivatives as computed by this routine.
            class(*), intent(inout), optional :: args
                !! An optional argument allowing the passing of data in/out of
                !! this routine.
        end subroutine
    end interface

    interface frequency_sweep
        module procedure :: frf_sweep_1
        module procedure :: frf_sweep_2
    end interface

    type frf_arg_container
        procedure(harmonic_ode), pointer, nopass :: fcn
        real(real64) :: frequency
        logical :: uses_optional_args
        class(*), allocatable :: optional_args
    end type

contains
! ------------------------------------------------------------------------------
    pure elemental function chirp(t, amp, span, f1Hz, f2Hz) result(rst)
        !! Evaluates a linear chirp function.
        !! The instantaneous frequency varies linearly,
        !! $$ f(t)=f_1+\frac{f_2-f_1}{T}t, $$
        !! giving phase \(\phi(t)=2\pi(f_1t+(f_2-f_1)t^2/(2T))\) and
        !! response \(x(t)=A\sin(\phi(t))\).
        real(real64), intent(in) :: t
            !! The value of the independent variable at which to evaluate the
            !! chirp.
        real(real64), intent(in) :: amp
            !! The amplitude.
        real(real64), intent(in) :: span
            !! The duration of the time it takes to sweep from the start
            !! frequency to the end frequency.
        real(real64), intent(in) :: f1Hz
            !! The lower excitation frequency, in Hz.
        real(real64), intent(in) :: f2Hz
            !! The upper excitation frequency, in Hz.
        real(real64) :: rst
            !! The value of the function at t.

        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
        real(real64) :: c
        c = (f2Hz - f1Hz) / span
        rst = amp * sin(2.0d0 * pi * t * (0.5d0 * c * t + f1Hz))
    end function
! ******************************************************************************
! HARMONIC_ODE_CONTAINER ROUTINES
! ------------------------------------------------------------------------------
    function frf_sweep_1(fcn, freq, iv, solver, ncycles, ntransient, &
        points, inHz, args) result(rst)
        !! Computes the frequency response of each equation of a system of
        !! harmonically excited ODE's by sweeping through frequency. 
        !!
        !! The amplitude and phase are determined by means of a harmonic 
        !! projection method.
        !!
        !! For uniformly sampled data, the retained complex response is
        !! $$ Y = b + j a, \qquad |Y| = \sqrt{a^{2} + b^{2}}, $$
        !! with phase $$ \phi = \operatorname{atan2}(a,b). $$
        !!
        !! where
        !! $$ a = \frac{2}{N} \sum_{k=0}^{N-1} y_k \cos(\omega t_k), \qquad
        !! b = \frac{2}{N} \sum_{k=0}^{N-1} y_k \sin(\omega t_k). $$
        use diffeq, only : runge_kutta_45
        use dynamics_error_handling
        procedure(harmonic_ode), pointer, intent(in) :: fcn
            !! A pointer to the routine containing the ODE's to integrate.
        real(real64), intent(in), dimension(:) :: freq
            !! An M-element array containing the frequency points at which the 
            !! solution should be computed.  Notice, whatever units are utilized
            !! for this array are also the units of the excitation_frequency
            !! property in fcn.  Additionally, this array cannot contain any
            !! zero-valued elements as the ODE solution time for each frequency 
            !! is determined by the period of oscillation and number of cycles.
        real(real64), intent(in), dimension(:) :: iv
            !! An N-element array containing the initial conditions for each of 
            !! the N ODEs.
        class(ode_integrator), intent(inout), optional, target :: solver
            !! An optional differential equation solver.  The default solver
            !! is the Dormand-Prince Runge-Kutta integrator from the DIFFEQ
            !! library.
        integer(int32), intent(in), optional :: ncycles
            !! An optional parameter controlling the number of cycles to 
            !! analyze when determining the amplitude and phase of the response.
            !! The default is 5.
        integer(int32), intent(in), optional :: ntransient
            !! An optional parameter controlling how many of the initial 
            !! "transient" cycles to ignore.  The default is 30.
        integer(int32), intent(in), optional :: points
            !! An optional parameter controlling how many evenly spaced 
            !! solution points should be considered per cycle.  The default is 
            !! 1000.
        logical, intent(in), optional :: inHz
            !! Set to true if the units of the frequency vector are in Hz.  If
            !! false, the units are assumed as rad/s.  The default is false such
            !! that the frequency units are assumed to be rad/s.
        class(*), intent(inout), optional :: args
            !! An optional argument allowing for passing of data in/out of the
            !! fcn subroutine.
        type(frf) :: rst
            !! The resulting frequency responses.

        ! Parameters
        real(real64), parameter :: zerotol = sqrt(epsilon(0.0d0))
        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

        ! Local Variables
        logical :: hz
        integer(int32) :: i, j, nfreq, neqn, nc, nt, ntotal, npts, ppc, &
            i1, ncpts
        real(real64) :: dt, omega, f
        real(real64), allocatable, dimension(:) :: ic, t
        real(real64), allocatable, dimension(:,:) :: sol
        type(ode_container) :: sys
        class(ode_integrator), pointer :: integrator
        type(runge_kutta_45), target :: default_integrator
        type(frf_arg_container) :: container
        
        ! Initialization
        if (present(ncycles)) then
            nc = ncycles
        else
            nc = 5
        end if
        if (present(ntransient)) then
            nt = ntransient
        else
            nt = 30
        end if
        if (present(points)) then
            ppc = points
        else
            ppc = 1000
        end if
        if (present(inHz)) then
            hz = inHz
        else
            hz = .false.
        end if
        nfreq = size(freq)
        neqn = size(iv)
        ntotal = nt + nc
        npts = ntotal * ppc
        ncpts = nc * ppc
        i1 = npts - ncpts + 1
        sys%fcn => sweep_eom

        ! Set up the optional argument container
        container%fcn => fcn
        container%uses_optional_args = present(args)
        if (present(args)) then
            allocate(container%optional_args, source = args)
        end if

        ! Set up the integrator
        if (present(solver)) then
            integrator => solver
        else
            integrator => default_integrator
        end if

        ! Input Checking
        if (.not.associated(fcn)) error stop DYN_NULL_POINTER_ERROR
        if (nfreq < 1 .or. neqn < 1) error stop DYN_INVALID_INPUT_ERROR
        if (nc < 1) error stop DYN_INVALID_INPUT_ERROR
        if (nt < 1) error stop DYN_INVALID_INPUT_ERROR
        if (ppc < 2) error stop DYN_INVALID_INPUT_ERROR
        do i = 1, nfreq
            if (abs(freq(i)) < zerotol) error stop DYN_ZERO_VALUED_FREQUENCY_ERROR
        end do

        ! Local Memory Allocation
        allocate(rst%responses(nfreq, neqn))
        allocate(rst%frequency(nfreq), source = freq)
        allocate(ic(neqn), source = iv)
        allocate(t(npts))

        ! Cycle over each frequency point
        do i = 1, nfreq
            ! Define the time vector
            if (hz) then
                f = freq(i)
                omega = 2.0d0 * pi * f
            else
                omega = freq(i)
                f = omega / (2.0d0 * pi)
            end if
            dt = (1.0d0 / f) / (ppc - 1.0d0)
            t = (/ (dt * j, j = 0, npts - 1) /)

            ! Set the frequency
            container%frequency = freq(i)

            ! Compute the solution
            call integrator%solve(sys, t, ic, args = container)
            sol = integrator%get_solution()

            ! Reset the initial conditions to the last solution point
            ic = sol(npts, 2:)

            ! Determine the magnitude and phase for each equation
            do j = 1, neqn
                rst%responses(i,j) = &
                    harmonic_projection(omega, sol(i1:npts,1), &
                    sol(i1:npts,j+1))
            end do

            ! Clear the solution buffer for the next time around
            call integrator%clear_buffer()
        end do
    end function

! ----------
    subroutine sweep_eom(x, y, dydx, args)
        real(real64), intent(in) :: x
        real(real64), intent(in), dimension(:) :: y
        real(real64), intent(out), dimension(:) :: dydx
        class(*), intent(inout), optional :: args

        select type (args)
        class is (frf_arg_container)
            if (args%uses_optional_args) then
                call args%fcn(args%frequency, x, y, dydx, args%optional_args)
            else
                call args%fcn(args%frequency, x, y, dydx)
            end if
        end select
    end subroutine

! ----------
    pure function harmonic_projection(omega, t, y) result(rst)
        !! Utilizes a harmonic projection method to determine the amplitude and
        !! phase of the supplied solution.
        real(real64), intent(in) :: omega
            !! The excitation frequency, in rad/s.
        real(real64), intent(in), dimension(:) :: t
            !! The solution time points.
        real(real64), intent(in), dimension(size(t)) :: y
            !! The solution vector.
        complex(real64) :: rst

        ! Local Variables
        integer(int32) :: i, n
        real(real64) :: a, b, r, phi

        ! Initialization
        n = size(t)
        a = 0.0d0
        b = 0.0d0

        ! Process
        do i = 1, n
            a = a + y(i) * cos(omega * t(i))
            b = b + y(i) * sin(omega * t(i))
        end do
        a = 2.0d0 * a / n
        b = 2.0d0 * b / n
        r = sqrt(a**2 + b**2)
        phi = atan2(a, b)
        rst = cmplx(r * cos(phi), r * sin(phi))
    end function

! ------------------------------------------------------------------------------
    function frf_sweep_2(fcn, nfreq, freq1, freq2, iv, solver, ncycles, &
        ntransient, points, inHz, args) result(rst)
        !! Computes the frequency response of each equation of a system of
        !! harmonically excited ODE's by sweeping through frequency.
        !!
        !! The amplitude and phase are determined by means of a harmonic 
        !! projection method.
        !!
        !! For uniformly sampled data, the retained complex response is
        !! $$ Y = b + j a, \qquad |Y| = \sqrt{a^{2} + b^{2}}, $$
        !! with phase $$ \phi = \operatorname{atan2}(a,b). $$
        !!
        !! where
        !! $$ a = \frac{2}{N} \sum_{k=0}^{N-1} y_k \cos(\omega t_k), \qquad
        !! b = \frac{2}{N} \sum_{k=0}^{N-1} y_k \sin(\omega t_k). $$
        procedure(harmonic_ode), pointer, intent(in) :: fcn
            !! A pointer to the routine containing the ODE's to integrate.
        integer(int32), intent(in) :: nfreq
            !! The number of frequency values to analyze.  This value must be
            !! at least 2.
        real(real64), intent(in) :: freq1
            !! The starting frequency.
        real(real64), intent(in) :: freq2
            !! The ending frequency.
        real(real64), intent(in), dimension(:) :: iv
            !! An N-element array containing the initial conditions for each of 
            !! the N ODEs.
        class(ode_integrator), intent(inout), optional, target :: solver
            !! An optional differential equation solver.  The default solver
            !! is the Dormand-Prince Runge-Kutta integrator from the DIFFEQ
            !! library.
        integer(int32), intent(in), optional :: ncycles
            !! An optional parameter controlling the number of cycles to 
            !! analyze when determining the amplitude and phase of the response.
            !! The default is 5.
        integer(int32), intent(in), optional :: ntransient
            !! An optional parameter controlling how many of the initial 
            !! "transient" cycles to ignore.  The default is 30.
        integer(int32), intent(in), optional :: points
            !! An optional parameter controlling how many evenly spaced 
            !! solution points should be considered per cycle.  The default is 
            !! 1000.
        logical, intent(in), optional :: inHz
            !! Set to true if the units of the frequency units are in Hz.  If
            !! false, the units are assumed as rad/s.  The default is false such
            !! that the frequency units are assumed to be rad/s.
        class(*), intent(inout), optional :: args
            !! An optional argument allowing for passing of data in/out of the
            !! fcn subroutine.
        type(frf) :: rst
            !! The resulting frequency responses.

        ! Local Variables
        integer(int32) :: i
        real(real64) :: df
        real(real64), allocatable, dimension(:) :: freq

        ! Input Checking
        if (.not.associated(fcn)) error stop DYN_NULL_POINTER_ERROR
        if (abs(freq1 - freq2) < sqrt(epsilon(freq1))) error stop DYN_INVALID_INPUT_ERROR
        if (nfreq < 2) error stop DYN_INVALID_INPUT_ERROR

        ! Process
        df = (freq2 - freq1) / (nfreq - 1.0d0)
        allocate(freq(nfreq))
        freq = (/ (df * i + freq1, i = 0, nfreq - 1) /)
        rst = frf_sweep_1(fcn, freq, iv, solver, ncycles, ntransient, &
            points, inHz = inHz, args = args)
    end function

end module
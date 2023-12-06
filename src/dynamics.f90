module dynamics
    use iso_fortran_env
    use diffeq
    use ferror
    implicit none
    private
    public :: DYN_MEMORY_ERROR
    public :: DYN_NULL_POINTER_ERROR
    public :: ode_excite
    public :: forced_ode_container
    public :: harmonic_ode_container
    public :: chirp
    public :: frequency_response

    integer(int32), parameter :: DYN_MEMORY_ERROR = DIFFEQ_MEMORY_ALLOCATION_ERROR
        !! Defines an error associated with memory allocations.
    integer(int32), parameter :: DYN_NULL_POINTER_ERROR = DIFFEQ_NULL_POINTER_ERROR
        !! Defines an error associated with a null pointer.
    integer(int32), parameter :: DYN_INVALID_INPUT_ERROR = DIFFEQ_INVALID_INPUT_ERROR
        !! Defines an error associated with an invalid input.

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
    end interface

    type, extends(ode_container) :: forced_ode_container
        !! An extension of the ode_container type from the DIFFEQ library
        !! that allows for the definition and analysis of forced systems of
        !! ODE's.
        procedure(ode_excite), pointer, nopass :: forcing_function
            !! A pointer to a routine defining the forcing function.
    end type

    type, extends(ode_container) :: harmonic_ode_container
        !! An extension of the ode_container type from the DIFFEQ library
        !! that allows for the definition and analysis of systems of ODE's
        !! exposed to harmonic excitation.
        real(real64) :: excitation_frequency
            !! The excitation frequency.
    contains
        procedure, public :: frequency_sweep => hoc_frf_sweep
    end type

    interface frequency_response
        !! Computes the frequency response functions for a system of ODE's.
        module procedure :: frf_fft
    end interface

contains
! ------------------------------------------------------------------------------
    pure elemental function chirp(t, amp, span, f1Hz, f2Hz) result(rst)
        !! Evaluates a linear chirp function.
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

! ------------------------------------------------------------------------------
    subroutine frf_fft(sys, span, iv, rst, fs, solver, win, freq, err)
        use spectrum
        !! Computes the frequency response of each equation in a system of 
        !! forced differential equations.
        class(forced_ode_container), intent(inout) :: sys
            !! The forced_ode_container object containing the equations to
            !! analyze.
        real(real64), intent(in) :: span
            !! The duration of the time-domain analysis.
        real(real64), intent(in), dimension(:) :: iv
            !! An N-element containing the initial conditions for each of the
            !! N differential equations being analyzed.
        complex(real64), intent(out), allocatable, dimension(:,:) :: rst
            !! An M-by-N matrix containing the response of each of the N 
            !! equations for M frequency values, where M is determined by the
            !! window size.
        real(real64), intent(in), optional :: fs
            !! An optional rate at which to sample the differential equation 
            !! solution.  The default rate is 1024 Hz, assuming that the input
            !! span is provided in units of seconds.
        class(ode_integrator), intent(inout), optional, target :: solver
            !! An optional differential equation solver.  The default solver
            !! is the Dormand-Prince Runge-Kutta integrator from the DIFFEQ
            !! library.
        class(window), intent(in), optional, target :: win
            !! An optional parameter allowing the differential equation solution
            !! to be windowed as part of the FFT process used to construct the
            !! system transfer functions.  The default is a Hamming window sized
            !! to contain the entirety of the solution.
        real(real64), intent(out), optional, allocatable, dimension(:) :: freq
            !! An optional array, that if provided, returns the M frequency 
            !! values, in the units of the input fs, where the frequency
            !! responses were calculated.
        class(errors), intent(inout), optional, target :: err
            !! An optional errors-based object that if provided 
            !! can be used to retrieve information relating to any errors 
            !! encountered during execution. If not provided, a default 
            !! implementation of the errors class is used internally to provide 
            !! error handling. Possible errors and warning messages that may be 
            !! encountered are as follows.
            !!
            !! - DYN_MEMORY_ERROR: Occurs if there are issues allocating memory.
            !! - DYN_NULL_POINTER_ERROR: Occurs if a null pointer is supplied.

        ! Local Variables
        integer(int32) :: i, npts, neqn, nfreq, flag
        real(real64) :: dt, sampleRate, df
        real(real64), allocatable, dimension(:) :: t, force
        real(real64), allocatable, dimension(:,:) :: sol
        class(ode_integrator), pointer :: integrator
        type(dprk45_integrator), target :: defaultIntegrator
        class(window), pointer :: w
        type(hamming_window), target :: defaultWindow
        class(errors), pointer :: errmgr
        type(errors), target :: deferr

        ! Initialization
        if (present(fs)) then
            sampleRate = fs
        else
            sampleRate = 1.024d3
        end if
        if (present(solver)) then
            integrator => solver
        else
            integrator => defaultIntegrator
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        dt = 1.0d0 / sampleRate
        npts = floor(span / dt) + 1
        neqn = size(iv)
        if (present(win)) then
            w => win
        else
            w => defaultWindow
            w%size = npts
        end if
        nfreq = compute_transform_length(w%size)

        ! Input Checking
        if (.not.associated(sys%forcing_function)) go to 10

        ! Memory Allocation
        allocate(rst(nfreq, neqn), stat = flag, source = (0.0d0, 0.0d0))
        if (flag == 0) allocate(t(npts), stat = flag)
        if (flag == 0) allocate(force(npts), stat = flag)
        if (flag == 0 .and. present(freq)) allocate(freq(nfreq), stat = flag)
        if (flag /= 0) go to 20

        ! Construct the array and forcing function
        do i = 1, npts
            t(i) = dt * (i - 1.0d0)
            force(i) = sys%forcing_function(t(i))
        end do

        ! Solve the ODE's
        sol = integrator%solve(sys, t, iv, errmgr)
        if (errmgr%has_error_occurred()) return

        ! Compute the transfer function of each result
        do i = 1, neqn
            rst(:,i) = siso_transfer_function(w, force, sol(:,i+1), &
                err = errmgr)
            if (errmgr%has_error_occurred()) return
        end do

        ! Build the frequency vector
        if (present(freq)) then
            df = frequency_bin_width(sampleRate, w%size)
            freq = (/ (df * i, i = 0, nfreq - 1) /)
        end if

        ! End
        return

        ! Null Forcing Function
    10  continue
        ! TO DO: Handle the error
        return

        ! Memory Allocation Error
    20  continue
        ! TO DO: Handle the error
        return
    end subroutine

! ------------------------------------------------------------------------------
    function hoc_frf_sweep(sys, freq, iv, solver, ncycles, ntransient, &
        points, err) result(rst)
        use spectrum, only : next_power_of_two
        !! Computes the frequency response of each equation of a system of
        !! harmonically excited ODE's by sweeping through frequency.
        class(harmonic_ode_container), intent(inout) :: sys
            !! The harmonic_ode_container object containing the equations to 
            !! analyze.  To properly use this object, extend the 
            !! harmonic_ode_container object and overload the ode routine to 
            !! define the ODE's.  Use the excitation_frequency property to
            !! obtain the desired frequency from the solver.
        real(real64), intent(in), dimension(:) :: freq
            !! An M-element array containing the frequency points at which the 
            !! solution should be computed.  Notice, whatever units are utilized
            !! for this array are also the units of the excitation_frequency
            !! property in @p sys.  It is recommended that the units be set to 
            !! Hz.  Additionally, this array cannot contain any zero-valued 
            !! elements as the ODE solution time for each frequency is 
            !! determined by the period of oscillation and number of cycles.
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
            !! The default is 20.
        integer(int32), intent(in), optional :: ntransient
            !! An optional parameter controlling how many of the initial 
            !! "transient" cycles to ignore.  The default is 200.
        integer(int32), intent(in), optional :: points
            !! An optional parameter controlling how many evenly spaced 
            !! solution points should be considered per cycle.  The default is 
            !! 1000.  Notice, there must be at least 2 points per cycle for the
            !! analysis to be effective.  The algorithm utilizes a discrete 
            !! Fourier transform to determine the phase and amplitude, and in 
            !! order to satisfy Nyquist conditions, the value must be at least 
            !! 2.
        class(errors), intent(inout), optional, target :: err
            !! An optional errors-based object that if provided 
            !! can be used to retrieve information relating to any errors 
            !! encountered during execution. If not provided, a default 
            !! implementation of the errors class is used internally to provide 
            !! error handling. Possible errors and warning messages that may be 
            !! encountered are as follows.
            !!
            !! - DYN_MEMORY_ERROR: Occurs if there are issues allocating memory.
            !! - DYN_NULL_POINTER_ERROR: Occurs if a null pointer is supplied.
            !! - DYN_INVALID_INPUT_ERROR: Occurs if an invalid parameter
            !!      is given.
        complex(real64), allocatable, dimension(:,:) :: rst
            !! The results matrix that will be allocated to M-by-N in size with
            !! each column containing the complex-valued frequency response of
            !! the corresponding ODE.

        ! Parameters
        real(real64), parameter :: zerotol = sqrt(epsilon(0.0d0))

        ! Local Variables
        integer(int32) :: i, j, nfreq, neqn, nc, nt, ntotal, npts, ppc, flag, &
            nfft, i1, ncpts
        real(real64) :: dt, tare, phase, amp
        real(real64), allocatable, dimension(:) :: ic, t
        real(real64), allocatable, dimension(:,:) :: sol
        complex(real64), allocatable, dimension(:) :: xpts
        class(ode_integrator), pointer :: integrator
        type(dprk45_integrator), target :: default_integrator
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        character(len = :), allocatable :: errmsg
        
        ! Initialization
        if (present(ncycles)) then
            nc = ncycles
        else
            nc = 20
        end if
        if (present(ntransient)) then
            nt = ntransient
        else
            nt = 200
        end if
        if (present(points)) then
            ppc = points
        else
            ppc = 1000
        end if
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        nfreq = size(freq)
        neqn = size(iv)
        ntotal = nt + nc
        npts = ntotal * ppc
        ncpts = nc * ppc
        i1 = npts - ncpts + 1
        nfft = 2**next_power_of_two(ppc * nc)

        ! Set up the integrator
        if (present(solver)) then
            integrator => solver
        else
            integrator => default_integrator
        end if

        ! Input Checking
        if (nc < 1) go to 20
        if (nt < 1) go to 30
        if (ppc < 2) go to 40
        do i = 1, nfreq
            if (abs(freq(i)) < zerotol) go to 50
        end do

        ! Local Memory Allocation
        allocate(rst(nfreq, neqn), stat = flag)
        if (flag == 0) allocate(ic(neqn), stat = flag, source = iv)
        if (flag == 0) allocate(t(npts), stat = flag)
        if (flag == 0) allocate(xpts(nfft), stat = flag, source = (0.0d0, 0.0d0))
        if (flag /= 0) go to 10

        ! Cycle over each frequency point
        do i = 1, nfreq
            ! Define the time vector
            dt = (1.0d0 / freq(i)) / (ppc - 1.0d0)
            t = (/ (dt * j, j = 0, npts - 1) /)

            ! Update the frequency
            sys%excitation_frequency = freq(i)

            ! Compute the solution
            sol = integrator%solve(sys, t, ic, errmgr)
            if (errmgr%has_error_occurred()) return

            ! Reset the initial conditions to the last solution point
            ic = sol(npts, 2:)

            ! Determine the magnitude and phase for each equation
            do j = 1, neqn
                rst(i,j) = get_magnitude_phase(sol(i1:,j+1), xpts)
            end do
        end do

        ! End
        return

        ! Memory Error
    10  continue
        allocate(character(len = 256) :: errmsg)
        write(errmsg, 100) "Memory allocation error flag ", flag, "."
        call err%report_error("frf_sweep", trim(errmsg), &
            DIFFEQ_MEMORY_ALLOCATION_ERROR)
        return

        ! Number of Cycles Error
    20  continue
        allocate(character(len = 256) :: errmsg)
        write(errmsg, 100) "The number of cycles to analyze must be at " // &
            "least 1; however, a value of ", nc, " was found."
        call errmgr%report_error("frf_sweep", trim(errmsg), &
            DYN_INVALID_INPUT_ERROR)
        return

        ! Number of Transient Cycles Error
    30  continue
        allocate(character(len = 256) :: errmsg)
        write(errmsg, 100) "The number of transient cycles must be at " // &
            "least 1; however, a value of ", nt, " was found."
        call errmgr%report_error("frf_sweep", trim(errmsg), &
            DYN_INVALID_INPUT_ERROR)
        return

        ! Points Per Cycle Error
    40  continue
        write(errmsg, 100) "The number of points per cycle must be at " // &
            "least 2; however, a value of ", ppc, " was found."
        call errmgr%report_error("frf_sweep", trim(errmsg), &
            DYN_INVALID_INPUT_ERROR)
        return

        ! Zero-Valued Frequency Error
    50  continue
        write(errmsg, 100) "A zero-valued frequency was found at index ", i, "."
        call errmgr%report_error("frf_sweep", trim(errmsg), &
            DYN_INVALID_INPUT_ERROR)
        return

        ! Formatting
    100 format(A, I0, A)
    end function

    ! ----------
    function get_magnitude_phase(x, xzeros) result(rst)
        use fftpack, only : fft
        use spectrum, only : compute_transform_length
        ! Arguments
        real(real64), intent(in), dimension(:) :: x
        complex(real64), intent(inout), dimension(:) :: xzeros
        complex(real64) :: rst

        ! Local Variables
        integer(int32) :: ind, m, n, nx
        real(real64) :: amp, phase

        ! Initialization
        nx = size(x)
        n = size(xzeros)
        m = compute_transform_length(n)
        amp = 0.5d0 * (maxval(x) - minval(x))

        ! Zero pad the data
        xzeros(:nx) = cmplx(x, 0.0d0, real64)
        xzeros(nx+1:) = cmplx(0.0d0, 0.0d0, real64)

        ! Compute the FFT to estimate the phase
        xzeros = fft(xzeros)
        ind = maxloc(abs(xzeros(:m)), 1)
        phase = atan2(aimag(xzeros(ind)), real(xzeros(ind)))
        rst = cmplx(amp * cos(phase), amp * sin(phase), real64)
    end function

! ------------------------------------------------------------------------------
end module
module dynamics
    use iso_fortran_env
    use diffeq
    use ferror
    use spectrum
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
        !! Computes the frequency response of each equation in a system of 
        !! forced differential equations.
        class(forced_ode_container), intent(inout) :: sys
            !! The forced_ode_container object.
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
            !!  can be used to retrieve information relating to any errors 
            !!  encountered during execution. If not provided, a default 
            !!  implementation of the errors class is used internally to provide 
            !!  error handling. Possible errors and warning messages that may be 
            !!  encountered are as follows.
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

! ------------------------------------------------------------------------------
end module
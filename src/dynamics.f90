module dynamics
    use iso_fortran_env
    use ferror
    use diffeq, only : ode_container, ode_integrator
    implicit none
    private
    public :: ode_excite
    public :: modal_excite
    public :: forced_ode_container
    public :: harmonic_ode_container
    public :: frf
    public :: chirp
    public :: frequency_response
    public :: compute_modal_damping
    public :: modal_response
    public :: normalize_mode_shapes
    public :: rotate_x

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

        subroutine modal_excite(freq, frc)
            !! Defines the interface to a routine for defining the forcing
            !! function for a modal frequency analysis.
            use iso_fortran_env, only : real64
            real(real64), intent(in) :: freq
                !! The excitation frequency.
            complex(real64), intent(out), dimension(:) :: frc
                !! An N-element array where the forcing function should be
                !! written.
        end subroutine
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
        generic, public :: frequency_sweep => hoc_frf_sweep, hoc_frf_sweep_2
        procedure, private :: hoc_frf_sweep
        procedure, private :: hoc_frf_sweep_2
    end type

    type frf
        real(real64), allocatable, dimension(:) :: frequency
        complex(real64), allocatable, dimension(:,:) :: responses
    end type

    interface frequency_response
        !! Computes the frequency response functions for a system of ODE's.
        module procedure :: frf_fft
        module procedure :: frf_modal_prop_damp
        module procedure :: frf_modal_prop_damp_2
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
    function frf_fft(sys, span, iv, fs, solver, win, err) result(rst)
        !! Computes the frequency response of each equation in a system of 
        !! forced differential equations.
        use spectrum
        use diffeq, only : dprk45_integrator
        use dynamics_error_handling
        class(forced_ode_container), intent(inout) :: sys
            !! The forced_ode_container object containing the equations to
            !! analyze.
        real(real64), intent(in) :: span
            !! The duration of the time-domain analysis.
        real(real64), intent(in), dimension(:) :: iv
            !! An N-element containing the initial conditions for each of the
            !! N differential equations being analyzed.
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
        type(frf) :: rst
            !! The resulting frequency responses.

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
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
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
        allocate(rst%responses(nfreq, neqn), stat = flag)
        if (flag == 0) allocate(rst%frequency(nfreq), stat = flag)
        if (flag == 0) allocate(t(npts), stat = flag)
        if (flag == 0) allocate(force(npts), stat = flag)
        if (flag /= 0) go to 20

        ! Construct the array and forcing function
        do i = 1, npts
            t(i) = dt * (i - 1.0d0)
            force(i) = sys%forcing_function(t(i))
        end do

        ! Solve the ODE's
        sol = integrator%solve(sys, t, iv, err = errmgr)
        if (errmgr%has_error_occurred()) return

        ! Compute the transfer function of each result
        do i = 1, neqn
            rst%responses(:,i) = siso_transfer_function(w, force, sol(:,i+1), &
                err = errmgr)
            if (errmgr%has_error_occurred()) return
        end do

        ! Build the frequency vector
        df = frequency_bin_width(sampleRate, w%size)
        rst%frequency = (/ (df * i, i = 0, nfreq - 1) /)

        ! End
        return

        ! Null Forcing Function
    10  continue
        call report_null_forcing_routine_error("frf_fft", errmgr)
        return

        ! Memory Allocation Error
    20  continue
        call report_memory_error("frf_fft", flag, errmgr)
        return
    end function

! ------------------------------------------------------------------------------
    function frf_modal_prop_damp(mass, stiff, alpha, beta, freq, frc, &
        modes, modeshapes, err) result(rst)
        !! Computes the frequency response functions for a 
        !! multi-degree-of-freedom system that uses proportional damping such
        !! that the damping matrix \( C \) is related to the stiffness an mass
        !! matrices by proportional damping coefficients \( \alpha \) and
        !! \( \beta \) by \( C = \alpha M + \beta K \).
        use linalg, only : eigen, sort, mtx_mult, LA_NO_OPERATION, LA_TRANSPOSE
        use dynamics_error_handling
        real(real64), intent(in), dimension(:,:) :: mass
            !! The N-by-N mass matrix for the system.  This matrix must be
            !! symmetric.
        real(real64), intent(in), dimension(:,:) :: stiff
            !! The N-by-N stiffness matrix for the system.  This matrix must be
            !! symmetric.
        real(real64), intent(in) :: alpha
            !! The mass damping factor, \( \alpha \).
        real(real64), intent(in) :: beta
            !! The stiffness damping factor, \( \beta \).
        real(real64), intent(in), dimension(:) :: freq
            !! An M-element array of frequency values at which to evaluate the
            !! frequency response functions.
        procedure(modal_excite), pointer, intent(in) :: frc
            !! A pointer to a routine used to compute the modal forcing 
            !! function.
        real(real64), intent(out), allocatable, optional, &
            dimension(:) :: modes
            !! An optional N-element allocatable array that, if supplied, will
            !! be used to retrieve the modal frequencies, in units of Hz.
        real(real64), intent(out), allocatable, optional, &
            dimension(:,:) :: modeshapes
            !! An optional N-by-N allocatable matrix that, if supplied, will be
            !! used to retrieve the N mode shapes with each vector occupying
            !! its own column.
        class(errors), intent(inout), optional, target :: err
            !! An optional errors-based object that if provided 
            !! can be used to retrieve information relating to any errors 
            !! encountered during execution. If not provided, a default 
            !! implementation of the errors class is used internally to provide 
            !! error handling. Possible errors and warning messages that may be 
            !! encountered are as follows.
            !!
            !! - DYN_MEMORY_ERROR: Occurs if there are issues allocating memory.
            !! - DYN_MATRIX_SIZE_ERROR: Occurs if the mass or stiffness matrices
            !!      are not square, or if the mass and stiffness matrices are
            !!      different sized.
            !! - DYN_NULL_POINTER_ERROR: Occurs if the forcing function pointer
            !!      is undefined.
        type(frf) :: rst
            !! The resulting frequency responses.

        ! Parameters
        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
        complex(real64), parameter :: j = (0.0d0, 1.0d0)
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        integer(int32) :: i, m, n, flag
        complex(real64) :: s
        real(real64), allocatable, dimension(:) :: lambda, zeta
        real(real64), allocatable, dimension(:,:) :: mmtx, kmtx
        complex(real64), allocatable, dimension(:) :: vals, q, f, u
        complex(real64), allocatable, dimension(:,:) :: vecs
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        m = size(freq)
        n = size(mass, 1)

        ! Input Checking
        if (size(mass, 2) /= n) go to 20
        if (size(stiff,1) /= size(stiff, 2)) go to 30
        if (size(stiff, 1) /= n .or. size(stiff, 2) /= n) go to 40
        if (.not.associated(frc)) go to 50

        ! TO DO: Check for symmetry

        ! Memory allocations
        allocate(mmtx(n, n), source = mass, stat = flag)
        if (flag == 0) allocate(kmtx(n, n), source = stiff, stat = flag)
        if (flag == 0) allocate(zeta(n), stat = flag)
        if (flag == 0) allocate(q(n), stat = flag)
        if (flag == 0) allocate(vals(n), stat = flag)
        if (flag == 0) allocate(f(n), stat = flag, source = zero)
        if (flag == 0) allocate(u(n), stat = flag)
        if (flag == 0) allocate(vecs(n, n), stat = flag)
        if (flag == 0) allocate(rst%responses(m, n), stat = flag)
        if (flag == 0) allocate(rst%frequency(m), source = freq, stat = flag)
        if (flag /= 0) go to 10

        ! Compute the eigenvalues and eigenvectors
        call eigen(kmtx, mmtx, vals, vecs = vecs, err = errmgr)
        if (errmgr%has_error_occurred()) return

        allocate(lambda(n), source = real(vals), stat = flag)
        if (flag /= 0) go to 10

        ! Compute the damping terms
        zeta = compute_modal_damping(lambda, alpha, beta)

        ! Compute each transfer function
        do i = 1, m
            call frc(freq(i), f)
            call mtx_mult(LA_TRANSPOSE, one, vecs, f, zero, u)
            s = j * (2.0d0 * pi * freq(i))
            q = u / (s**2 + 2.0d0 * zeta * sqrt(lambda) * s + lambda)
            call mtx_mult(LA_NO_OPERATION, one, vecs, q, zero, rst%responses(i,:))
        end do

        ! If needed, return the modal frequencies and mode shapes
        if (present(modes) .or. present(modeshapes)) then
            ! Sort the modal information
            call sort(vals, vecs)
        end if

        if (present(modes)) then
            allocate(modes(n), source = sqrt(real(vals)) / (2.0d0 * pi), &
                stat = flag)
            if (flag /= 0) go to 10
        end if

        if (present(modeshapes)) then
            allocate(modeshapes(n, n), source = real(vecs), stat = flag)
            if (flag /= 0) go to 10
        end if

        ! End
        return

        ! Memory error
    10  continue
        call report_memory_error("frf_modal_prop_damp", flag, errmgr)
        return

        ! Error: Mass matrix is not square
    20  continue
        call report_nonsquare_mass_matrix_error("frf_modal_prop_damp", &
            size(mass, 1), size(mass, 2), errmgr)
        return

        ! Error: Stiffness matrix is not square
    30  continue
        call report_nonsquare_stiffness_matrix_error("frf_modal_prop_damp", &
            size(stiff, 1), size(stiff, 2), errmgr)
        return

        ! Error: Stiffness matrix is not sized correctly
    40  continue
        call report_matrix_size_mismatch_error("frf_modal_prop_damp", &
            "mass", "stiffness", size(mass, 1), size(mass, 2), &
            size(stiff, 1), size(stiff, 2), errmgr)
        return

        ! Null forcing term pointer
    50  continue
        call report_null_forcing_routine_error("frf_modal_prop_damp", &
            errmgr)
        return
    end function

! ------------------------------------------------------------------------------
    function frf_modal_prop_damp_2(mass, stiff, alpha, beta, nfreq, freq1, &
        freq2, frc, modes, modeshapes, err) result(rst)
        !! Computes the frequency response functions for a 
        !! multi-degree-of-freedom system that uses proportional damping such
        !! that the damping matrix \( C \) is related to the stiffness an mass
        !! matrices by proportional damping coefficients \( \alpha \) and
        !! \( \beta \) by \( C = \alpha M + \beta K \).
        use linalg, only : eigen, sort, mtx_mult, LA_NO_OPERATION, LA_TRANSPOSE
        use dynamics_error_handling
        real(real64), intent(in), dimension(:,:) :: mass
            !! The N-by-N mass matrix for the system.  This matrix must be
            !! symmetric.
        real(real64), intent(in), dimension(:,:) :: stiff
            !! The N-by-N stiffness matrix for the system.  This matrix must be
            !! symmetric.
        real(real64), intent(in) :: alpha
            !! The mass damping factor, \( \alpha \).
        real(real64), intent(in) :: beta
            !! The stiffness damping factor, \( \beta \).
        integer(int32), intent(in) :: nfreq
            !! The number of frequency values to analyze.  This value must be
            !! at least 2.
        real(real64), intent(in) :: freq1
            !! The starting frequency.  It is recommended that the units be set
            !! to Hz.
        real(real64), intent(in) :: freq2
            !! The ending frequency.  It is recommended that the units be set to
            !! Hz.
        procedure(modal_excite), pointer, intent(in) :: frc
            !! A pointer to a routine used to compute the modal forcing 
            !! function.
        real(real64), intent(out), allocatable, optional, &
            dimension(:) :: modes
            !! An optional N-element allocatable array that, if supplied, will
            !! be used to retrieve the modal frequencies, in units of Hz.
        real(real64), intent(out), allocatable, optional, &
            dimension(:,:) :: modeshapes
            !! An optional N-by-N allocatable matrix that, if supplied, will be
            !! used to retrieve the N mode shapes with each vector occupying
            !! its own column.
        class(errors), intent(inout), optional, target :: err
            !! An optional errors-based object that if provided 
            !! can be used to retrieve information relating to any errors 
            !! encountered during execution. If not provided, a default 
            !! implementation of the errors class is used internally to provide 
            !! error handling. Possible errors and warning messages that may be 
            !! encountered are as follows.
            !!
            !! - DYN_MEMORY_ERROR: Occurs if there are issues allocating memory.
            !! - DYN_MATRIX_SIZE_ERROR: Occurs if the mass or stiffness matrices
            !!      are not square, or if the mass and stiffness matrices are
            !!      different sized.
            !! - DYN_NULL_POINTER_ERROR: Occurs if the forcing function pointer
            !!      is undefined.
        type(frf) :: rst
            !! The resulting frequency responses.

        ! Local Variables
        integer(int32) :: i, flag
        real(real64) :: df
        real(real64), allocatable, dimension(:) :: freq
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Checking
        if (abs(freq1 - freq2) < sqrt(epsilon(freq1))) then
            call report_zero_difference_error("frf_modal_prop_damp_2", &
                "freq1", freq1, "freq2", freq2, DYN_INVALID_INPUT_ERROR, &
                errmgr)
            return
        end if
        if (nfreq < 2) then
            call report_generic_counting_error("frf_modal_prop_damp_2", &
                "The number of frequency points must be at least 2, " // &
                "but was found to be ", nfreq, ".", DYN_INVALID_INPUT_ERROR, &
                errmgr)
            return
        end if

        ! Process
        df = (freq2 - freq1) / (nfreq - 1.0d0)
        allocate(freq(nfreq), stat = flag)
        if (flag /= 0) then
            call report_memory_error("frf_modal_prop_damp_2", flag, errmgr)
            return
        end if
        freq = (/ (df * i + freq1, i = 0, nfreq - 1) /)
        rst = frequency_response(mass, stiff, alpha, beta, freq, frc, modes, &
            modeshapes, err)
    end function

! ------------------------------------------------------------------------------
    pure elemental function compute_modal_damping(lambda, alpha, beta) &
        result(rst)
        !! Computes the modal damping factors (\ \zeta_i \) given the
        !! proportional damping terms (\ alpha \) and (\ beta \) where
        !! (\ \alpha + \beta \omega_{i}^2 = 2 \zeta_{i} \omega_{i} \) and
        !! (\ \lambda_{i} = \omega_{i}^2 \).
        real(real64), intent(in) :: lambda
            !! The square of the modal frequency - the eigen value.
        real(real64), intent(in) :: alpha
            !! The mass damping factor, \( \alpha \).
        real(real64), intent(in) :: beta
            !! The stiffness damping factor, \( \beta \).
        real(real64) rst
            !! The modal damping parameter.

        ! Local Variables
        integer(int32) :: n

        ! Process
        rst = (alpha + beta * lambda) / (2.0d0 * sqrt(lambda))
    end function

! ------------------------------------------------------------------------------
    subroutine modal_response(mass, stiff, freqs, modeshapes, err)
        !! Computes the modal frequencies and modes shapes for 
        !! multi-degree-of-freedom system.
        use dynamics_error_handling
        use linalg, only : eigen, sort
        real(real64), intent(in), dimension(:,:) :: mass
            !! The N-by-N mass matrix for the system.  This matrix must be
            !! symmetric.
        real(real64), intent(in), dimension(:,:) :: stiff
            !! The N-by-N stiffness matrix for the system.  This matrix must
            !! be symmetric.
        real(real64), intent(out), allocatable, dimension(:) :: freqs
            !! An allocatable N-element array where the modal frequencies will
            !! be returned in ascending order with units of Hz.
        real(real64), intent(out), allocatable, optional, dimension(:,:) :: &
            modeshapes
            !! An optional, allocatable N-by-N matrix where the N mode shapes
            !! for the system will be returned.  The mode shapes are stored in
            !! columns.
        class(errors), intent(inout), optional, target :: err

        ! Parameters
        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

        ! Local Variables
        integer(int32) :: n, flag
        real(real64), allocatable, dimension(:,:) :: mmtx, kmtx
        complex(real64), allocatable, dimension(:) :: vals
        complex(real64), allocatable, dimension(:,:) :: vecs
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
        n = size(mass, 1)

        ! Input Checking
        if (size(mass, 2) /= n) go to 10
        if (size(stiff, 1) /= size(stiff, 2)) go to 20
        if (size(stiff, 1) /= n .or. size(stiff, 2) /= n) go to 30

        ! TO DO: Check for symmetry

        ! Memory allocations
        allocate(mmtx(n, n), source = mass, stat = flag)
        if (flag == 0) allocate(kmtx(n, n), source = stiff, stat = flag)
        if (flag == 0) allocate(vals(n), stat = flag)
        if (flag == 0 .and. present(modeshapes)) &
            allocate(vecs(n, n), stat = flag)
        if (flag /= 0) go to 40

        ! Solve the eigen problem
        if (present(modeshapes)) then
            call eigen(kmtx, mmtx, vals, vecs = vecs, err = errmgr)
            if (errmgr%has_error_occurred()) return

            call sort(vals, vecs)
            allocate(modeshapes(n, n), source = real(vecs), stat = flag)
            if (flag /= 0) go to 40
        else
            call eigen(kmtx, mmtx, vals, err = errmgr)
            if (errmgr%has_error_occurred()) return
            call sort(vals)
        end if

        ! Convert the eigenvalues to frequency values
        allocate(freqs(n), source = sqrt(abs(real(vals))) / (2.0d0 * pi), &
            stat = flag)
        if (flag /= 0) go to 40

        ! End
        return

        ! Non-square mass matrix error handler
    10  continue
        call report_nonsquare_mass_matrix_error("modal_response", &
            size(mass, 1), size(mass, 2), errmgr)
        return

        ! Non-square stiffness matrix error handler
    20  continue
        call report_nonsquare_stiffness_matrix_error("modal_response", &
            size(stiff, 1), size(stiff, 2), errmgr)
        return

        ! Stiffness and mass matrix size mismatch error handler
    30  continue
        call report_matrix_size_mismatch_error("modal_response", &
            "mass", "stiffness", size(mass, 1), size(mass, 2), &
            size(stiff, 1), size(stiff, 2), errmgr)
        return

        ! Memory error handler
    40  continue
        call report_memory_error("modal_response", flag, errmgr)
        return
    end subroutine

! ------------------------------------------------------------------------------
    subroutine normalize_mode_shapes(x)
        !! Normalizes mode shape vectors such that the largest magnitude
        !! value in the vector is one.
        real(real64), intent(inout), dimension(:,:) :: x
            !! The matrix of mode shape vectors with one vector per column.

        ! Local Variables
        integer(int32) :: i, loc
        real(real64) :: factor

        ! Process
        do i = 1, size(x, 2)
            loc = maxloc(abs(x(:,i)), 1)
            factor = x(loc, i)
            x(:,i) = x(:,i) / factor
        end do
    end subroutine



! ******************************************************************************
! HARMONIC_ODE_CONTAINER ROUTINES
! ------------------------------------------------------------------------------
    function hoc_frf_sweep(sys, freq, iv, solver, ncycles, ntransient, &
        points, err) result(rst)
        !! Computes the frequency response of each equation of a system of
        !! harmonically excited ODE's by sweeping through frequency.
        use spectrum, only : next_power_of_two
        use diffeq, only : dprk45_integrator
        use dynamics_error_handling
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
            !! - DYN_ZERO_VALUED_FREQUENCY_ERROR: Occurs if a zero-valued 
            !!      frequency was supplied.
        type(frf) :: rst
            !! The resulting frequency responses.

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
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if
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
        allocate(rst%responses(nfreq, neqn), stat = flag)
        if (flag == 0) allocate(rst%frequency(nfreq), source = freq, &
            stat = flag)
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
            sol = integrator%solve(sys, t, ic, err = errmgr)
            if (errmgr%has_error_occurred()) return

            ! Reset the initial conditions to the last solution point
            ic = sol(npts, 2:)

            ! Determine the magnitude and phase for each equation
            do j = 1, neqn
                rst%responses(i,j) = get_magnitude_phase(sol(i1:,j+1), xpts)
            end do
        end do

        ! End
        return

        ! Memory Error
    10  continue
        call report_memory_error("hoc_frf_sweep", flag, errmgr)
        return

        ! Number of Cycles Error
    20  continue
        call report_generic_counting_error("hoc_frf_sweep", &
            "The number of cycles to analyze must be at least 1; " // &
            "however, a value of ", nc, " was found.", &
            DYN_INVALID_INPUT_ERROR, errmgr)
        return

        ! Number of Transient Cycles Error
    30  continue
        call report_generic_counting_error("hoc_frf_sweep", &
            "The number of transient cycles must be at least 1; " // &
            "however, a value of ", nt, " was found.", &
            DYN_INVALID_INPUT_ERROR, errmgr)
        return

        ! Points Per Cycle Error
    40  continue
        call report_generic_counting_error("hoc_frf_sweep", &
            "The number of points per cycle must be at least 2; " // &
            "however, a value of ", ppc, " was found.", &
            DYN_INVALID_INPUT_ERROR, errmgr)
        return

        ! Zero-Valued Frequency Error
    50  continue
        call report_zero_valued_frequency_error("hoc_frf_sweep", i, errmgr)
        return
    end function

    ! ----------
    function get_magnitude_phase(x, xzeros) result(rst)
        !! Returns the magnitude and phase of a signal.
        use fftpack, only : fft
        use spectrum, only : compute_transform_length
        ! Arguments
        real(real64), intent(in), dimension(:) :: x
            !! The array containing the signal.
        complex(real64), intent(inout), dimension(:) :: xzeros
            !! A workspace array for the FFT operation.
        complex(real64) :: rst
            !! The complex-valued result defining both magnitude and phase.

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
    function hoc_frf_sweep_2(sys, nfreq, freq1, freq2, iv, solver, ncycles, &
        ntransient, points, err) result(rst)
        !! Computes the frequency response of each equation of a system of
        !! harmonically excited ODE's by sweeping through frequency.
        use spectrum, only : next_power_of_two
        use diffeq, only : dprk45_integrator
        use dynamics_error_handling
        class(harmonic_ode_container), intent(inout) :: sys
            !! The harmonic_ode_container object containing the equations to 
            !! analyze.  To properly use this object, extend the 
            !! harmonic_ode_container object and overload the ode routine to 
            !! define the ODE's.  Use the excitation_frequency property to
            !! obtain the desired frequency from the solver.
        integer(int32), intent(in) :: nfreq
            !! The number of frequency values to analyze.  This value must be
            !! at least 2.
        real(real64), intent(in) :: freq1
            !! The starting frequency.  It is recommended that the units be set
            !! to Hz.
        real(real64), intent(in) :: freq2
            !! The ending frequency.  It is recommended that the units be set to
            !! Hz.
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
            !! - DYN_ZERO_VALUED_FREQUENCY_ERROR: Occurs if a zero-valued 
            !!      frequency was supplied.
        type(frf) :: rst
            !! The resulting frequency responses.

        ! Local Variables
        integer(int32) :: i, flag
        real(real64) :: df
        real(real64), allocatable, dimension(:) :: freq
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        
        ! Initialization
        if (present(err)) then
            errmgr => err
        else
            errmgr => deferr
        end if

        ! Input Checking
        if (abs(freq1 - freq2) < sqrt(epsilon(freq1))) then
            call report_zero_difference_error("hoc_frf_sweep_2", &
                "freq1", freq1, "freq2", freq2, DYN_INVALID_INPUT_ERROR, &
                errmgr)
            return
        end if
        if (nfreq < 2) then
            call report_generic_counting_error("hoc_frf_sweep_2", &
                "The number of frequency points must be at least 2, " // &
                "but was found to be ", nfreq, ".", DYN_INVALID_INPUT_ERROR, &
                errmgr)
            return
        end if

        ! Process
        df = (freq2 - freq1) / (nfreq - 1.0d0)
        allocate(freq(nfreq), stat = flag)
        if (flag /= 0) then
            call report_memory_error("hoc_frf_sweep_2", flag, errmgr)
            return
        end if
        freq = (/ (df * i + freq1, i = 0, nfreq - 1) /)
        rst = sys%frequency_sweep(freq, iv, solver, ncycles, ntransient, &
            points, err)
    end function



! ******************************************************************************
! KINEMATICS ROUTINES
! ******************************************************************************
    pure function rotate_x(angle) result(rst)
        !! Constructs the rotation matrix describing the rotation about an
        !! x-axis such that 
        !! \( \overrightarrow{r_2} = \textbf{R} \overrightarrow{r_1} \).
        !!
        !! $$ \textbf{R} = \begin{matrix} 1 & 0 & 0 \\ 0 & \cos{\theta_x} &
        !! -\sin{\theta_x} \\ 0 & \sin{\theta_x} & \cos{\theta_x} \\
        !! \end{matrix}
        real(real64), intent(in) :: angle
            !! The rotation angle, in radians.
        real(real64) :: rst(3, 3)
            !! The resulting 3-by-3 matrix.

        ! Local Variables
        real(real64) :: c, s

        ! Process
        c = cos(angle)
        s = sin(angle)
        rst = reshape([1.0d0, 0.0d0, 0.0d0, 0.0d0, c, s, 0.0d0, -s, c], [3, 3])
    end function

! ------------------------------------------------------------------------------
end module
module dynamics_frequency_response
    use iso_fortran_env
    use ferror
    use diffeq, only : ode_container, ode_integrator
    use dynamics_error_handling
    use spectrum
    use fstats
    implicit none
    private
    public :: ode_excite
    public :: modal_excite
    public :: harmonic_ode
    public :: frf
    public :: mimo_frf
    public :: chirp
    public :: frequency_response
    public :: frequency_sweep
    public :: compute_modal_damping
    public :: modal_response
    public :: normalize_mode_shapes
    public :: evaluate_accelerance_frf_model
    public :: evaluate_receptance_frf_model
    public :: fit_frf
    public :: FRF_ACCELERANCE_MODEL
    public :: FRF_RECEPTANCE_MODEL
    public :: regression_statistics
    public :: iteration_controls
    public :: lm_solver_options
    public :: convergence_info
    public :: ode_integrator

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
                !! The excitation frequency.  When used as a part of a frequency
                !! response calculation, this value will have the same units as
                !! the frequency values provided to the frequency response
                !! routine.
            complex(real64), intent(out), dimension(:) :: frc
                !! An N-element array where the forcing function should be
                !! written.
        end subroutine

        pure subroutine harmonic_ode(freq, t, x, dxdt, args)
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

    type frf
        !! A container for a frequency response function, or series of frequency
        !! response functions.
        real(real64), allocatable, dimension(:) :: frequency
            !! An N-element array containing the frequency values at which the 
            !! FRF is provided.  The units of this array are the same as the
            !! units of the frequency values passed to the routine used to 
            !! compute the frequency response.
        complex(real64), allocatable, dimension(:,:) :: responses
            !! An N-by-M matrix containing the M frequency response functions
            !! evaluated at each of the N frequency points.
    end type

    type mimo_frf
        !! A container for the frequency responses of a system of multiple 
        !! inputs and multiple outputs (MIMO).
        real(real64), allocatable, dimension(:) :: frequency
            !! An N-element array containing the frequency values at which the 
            !! FRF is provided.  The units of this array are the same as the
            !! units of the frequency values passed to the routine used to 
            !! compute the frequency response.
        complex(real64), allocatable, dimension(:,:,:) :: responses
            !! An N-by-M-by-P array containing the frequency response functions
            !! for each of the M outputs corresponding to each of the P inputs.
    end type

    interface frequency_response
        !! Computes the frequency response functions for a system of ODE's.
        module procedure :: frf_modal_prop_damp
        module procedure :: frf_modal_prop_damp_2
        module procedure :: siso_freqres
        module procedure :: mimo_freqres
    end interface

    interface frequency_sweep
        module procedure :: frf_sweep_1
        module procedure :: frf_sweep_2
    end interface

    interface evaluate_accelerance_frf_model
        module procedure :: evaluate_accelerance_frf_model_scalar
        module procedure :: evaluate_accelerance_frf_model_array
    end interface

    interface evaluate_receptance_frf_model
        module procedure :: evaluate_receptance_frf_model_scalar
        module procedure :: evaluate_receptance_frf_model_array
    end interface

! ------------------------------------------------------------------------------
    integer(int32), parameter :: FRF_ACCELERANCE_MODEL = 1
        !! Defines an accelerance frequency response model.
    integer(int32), parameter :: FRF_RECEPTANCE_MODEL = 2
        !! Defines a receptance frequency response model.

! ------------------------------------------------------------------------------
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
            !! frequency response functions, in units of rad/s.
        procedure(modal_excite), pointer, intent(in) :: frc
            !! A pointer to a routine used to compute the modal forcing 
            !! function.
        real(real64), intent(out), allocatable, optional, &
            dimension(:) :: modes
            !! An optional N-element allocatable array that, if supplied, will
            !! be used to retrieve the modal frequencies, in units of rad/s.
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
            s = j * freq(i)
            q = u / (s**2 + 2.0d0 * zeta * sqrt(lambda) * s + lambda)
            call mtx_mult(LA_NO_OPERATION, one, vecs, q, zero, rst%responses(i,:))
        end do

        ! If needed, return the modal frequencies and mode shapes
        if (present(modes) .or. present(modeshapes)) then
            ! Sort the modal information
            call sort(vals, vecs)
        end if

        if (present(modes)) then
            allocate(modes(n), source = sqrt(real(vals)), &
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
            !! The starting frequency, in units of rad/s.
        real(real64), intent(in) :: freq2
            !! The ending frequency, in units of rad/s.
        procedure(modal_excite), pointer, intent(in) :: frc
            !! A pointer to a routine used to compute the modal forcing 
            !! function.
        real(real64), intent(out), allocatable, optional, &
            dimension(:) :: modes
            !! An optional N-element allocatable array that, if supplied, will
            !! be used to retrieve the modal frequencies, in units of rad/s.
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
        !! Computes the modal damping factors \( \zeta_i \) given the
        !! proportional damping terms \( \alpha \) and \( \beta \) where
        !! \( \alpha + \beta \omega_{i}^2 = 2 \zeta_{i} \omega_{i} \),
        !! \( \lambda_{i} = \omega_{i}^2 \), and \( \lambda_i \) is the
        !! \( i^{th} \) eigenvalue of the system.
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
            !! be returned in ascending order with units of rad/s.
        real(real64), intent(out), allocatable, optional, dimension(:,:) :: &
            modeshapes
            !! An optional, allocatable N-by-N matrix where the N mode shapes
            !! for the system will be returned.  The mode shapes are stored in
            !! columns.
        class(errors), intent(inout), optional, target :: err

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
        allocate(freqs(n), source = sqrt(abs(real(vals))), &
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
    function frf_sweep_1(fcn, freq, iv, solver, ncycles, ntransient, &
        points, args, err) result(rst)
        !! Computes the frequency response of each equation of a system of
        !! harmonically excited ODE's by sweeping through frequency.
        use spectrum, only : next_power_of_two
        use fftpack, only : zffti
        use diffeq, only : runge_kutta_45
        use dynamics_error_handling
        procedure(harmonic_ode), pointer, intent(in) :: fcn
            !! A pointer to the routine containing the ODE's to integrate.
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
        class(*), intent(inout), optional :: args
            !! An optional argument allowing for passing of data in/out of the
            !! fcn subroutine.
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
            nfft, i1, ncpts, lsave
        real(real64) :: dt, tare, phase, amp
        real(real64), allocatable, dimension(:) :: ic, t, wsave
        real(real64), allocatable, dimension(:,:) :: sol
        complex(real64), allocatable, dimension(:) :: xpts
        type(ode_container) :: sys
        class(ode_integrator), pointer :: integrator
        type(runge_kutta_45), target :: default_integrator
        class(errors), pointer :: errmgr
        type(errors), target :: deferr
        type(frf_arg_container) :: container
        
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
        lsave = 4 * nfft + 15
        sys%fcn => sweep_eom

        ! Set up the optional argument container
        container%fcn => fcn
        container%uses_optional_args = present(args)
        if (present(args)) then
            container%optional_args = args
        end if

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
        if (flag == 0) allocate(wsave(lsave), stat = flag)
        if (flag /= 0) go to 10

        ! Set up the FFT calculations for the get_magnitude_phase routine
        call zffti(nfft, wsave)

        ! Cycle over each frequency point
        do i = 1, nfreq
            ! Define the time vector
            dt = (1.0d0 / freq(i)) / (ppc - 1.0d0)
            t = (/ (dt * j, j = 0, npts - 1) /)

            ! Set the frequency
            container%frequency = freq(i)

            ! Compute the solution
            call integrator%solve(sys, t, ic, args = container, err = errmgr)
            if (errmgr%has_error_occurred()) return
            sol = integrator%get_solution()

            ! Reset the initial conditions to the last solution point
            ic = sol(npts, 2:)

            ! Determine the magnitude and phase for each equation
            do j = 1, neqn
                rst%responses(i,j) = &
                    get_magnitude_phase(sol(i1:npts,j+1), xpts, wsave)
            end do

            ! Clear the solution buffer for the next time around
            call integrator%clear_buffer()
        end do

        ! End
        return

        ! Memory Error
    10  continue
        call report_memory_error("frf_sweep_1", flag, errmgr)
        return

        ! Number of Cycles Error
    20  continue
        call report_generic_counting_error("frf_sweep_1", &
            "The number of cycles to analyze must be at least 1; " // &
            "however, a value of ", nc, " was found.", &
            DYN_INVALID_INPUT_ERROR, errmgr)
        return

        ! Number of Transient Cycles Error
    30  continue
        call report_generic_counting_error("frf_sweep_1", &
            "The number of transient cycles must be at least 1; " // &
            "however, a value of ", nt, " was found.", &
            DYN_INVALID_INPUT_ERROR, errmgr)
        return

        ! Points Per Cycle Error
    40  continue
        call report_generic_counting_error("frf_sweep_1", &
            "The number of points per cycle must be at least 2; " // &
            "however, a value of ", ppc, " was found.", &
            DYN_INVALID_INPUT_ERROR, errmgr)
        return

        ! Zero-Valued Frequency Error
    50  continue
        call report_zero_valued_frequency_error("frf_sweep_1", i, errmgr)
        return
    end function

    ! ----------
    pure subroutine sweep_eom(x, y, dydx, args)
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
    function get_magnitude_phase(x, xzeros, wsave) result(rst)
        !! Returns the magnitude and phase of a signal.
        use fftpack, only : zfftf
        use spectrum, only : compute_transform_length
        ! Arguments
        real(real64), intent(in), dimension(:) :: x
            !! The array containing the signal.
        complex(real64), intent(out), dimension(:) :: xzeros
            !! A workspace array for the FFT operation.
        real(real64), intent(in), dimension(:) :: wsave
            !! A workspace array for the FFT operation
        complex(real64) :: rst
            !! The complex-valued result defining both magnitude and phase.

        ! Parameters
        complex(real64), parameter :: czero = (0.0d0, 0.0d0)

        ! Local Variables
        integer(int32) :: i, ind, m, n, nx

        ! Initialization
        nx = size(x)
        n = size(xzeros)
        m = compute_transform_length(n)

        ! Zero pad the data
        xzeros(:nx) = cmplx(x, 0.0d0, real64)
        xzeros(nx+1:) = czero

        ! Compute the FFT to estimate the phase
        call zfftf(n, xzeros, wsave)
        ind = maxloc(abs(xzeros(2:m)), 1) + 1 ! start at 2 to avoid any DC issues
        rst = 2.0d0 * xzeros(ind) / nx
    end function

! ------------------------------------------------------------------------------
    function frf_sweep_2(fcn, nfreq, freq1, freq2, iv, solver, ncycles, &
        ntransient, points, args, err) result(rst)
        !! Computes the frequency response of each equation of a system of
        !! harmonically excited ODE's by sweeping through frequency.
        use spectrum, only : next_power_of_two
        use diffeq, only : runge_kutta_45
        use dynamics_error_handling
        procedure(harmonic_ode), pointer, intent(in) :: fcn
            !! A pointer to the routine containing the ODE's to integrate.
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
        class(*), intent(inout), optional :: args
            !! An optional argument allowing for passing of data in/out of the
            !! fcn subroutine.
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
        rst = frf_sweep_1(fcn, freq, iv, solver, ncycles, ntransient, &
            points, args = args, err = err)
    end function

! ******************************************************************************
! VERSION 1.0.5 ADDITIONS
! ------------------------------------------------------------------------------
function siso_freqres(x, y, fs, win, method, err) result(rst)
    !! Estimates the frequency response of a single-input, single-output (SISO)
    !! system.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the excitation signal.
    real(real64), intent(in), dimension(:) :: y
        !! An N-element array containing the response signal.
    real(real64), intent(in) :: fs
        !! The sampling frequency, in Hz.
    class(window), intent(in), optional, target :: win
        !! The window to apply to the data.  If nothing is supplied, no window
        !! is applied.
    integer(int32), intent(in), optional :: method
        !! Enter 1 to utilize an H1 estimator; else, enter 2 to utilize an
        !! H2 estimator.  The default is an H1 estimator.
        !!
        !! An H1 estimator is defined as the cross-spectrum of the input and
        !! response signals divided by the energy spectral density of the input.
        !! An H2 estimator is defined as the energy spectral density of the
        !! response divided by the cross-spectrum of the input and response
        !! signals.
        !!
        !! $$ H_{1} = \frac{P_{xy}}{P_{xx}} $$
        !!
        !! $$ H_{2} = \frac{P_{yy}}{P_{xy}} $$
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to provide 
        !! error handling. Possible errors and warning messages that may be 
        !! encountered are as follows.
        !!
        !! - DYN_MEMORY_ERROR: Occurs if there are issues allocating memory.
        !!
        !! - DYN_ARRAY_SIZE_ERROR: Occurs if x and y are not the same size.
    type(frf) :: rst
        !! The resulting frequency response function.

    ! Local Variables
    integer(int32) :: i, npts, nfreq, meth, flag
    real(real64) :: df
    class(window), pointer :: wptr
    type(rectangular_window), target :: defwin
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(x)
    if (present(win)) then
        wptr => win
    else
        defwin%size = npts
        wptr => defwin
    end if
    if (present(method)) then
        if (method == 2) then
            meth = SPCTRM_H2_ESTIMATOR
        else
            meth = SPCTRM_H1_ESTIMATOR
        end if
    else
        meth = SPCTRM_H1_ESTIMATOR
    end if
    nfreq = compute_transform_length(wptr%size)
    allocate(rst%frequency(nfreq), stat = flag)
    if (flag == 0) allocate(rst%responses(nfreq, 1), stat = flag)
    if (flag /= 0) then
        call report_memory_error("siso_freqres", flag, errmgr)
        return
    end if

    ! Input Checking
    if (size(y) /= npts) then
        call report_array_size_error("siso_freqres", "y", npts, size(y), errmgr)
        return
    end if

    ! Compute the transfer function
    rst%responses(:,1) = siso_transfer_function(wptr, x, y, etype = meth, &
        err = errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the frequency vector
    df = frequency_bin_width(fs, wptr%size)
    rst%frequency = (/ (df * i, i = 0, nfreq - 1) /)
end function

! ------------------------------------------------------------------------------
function mimo_freqres(x, y, fs, win, method, err) result(rst)
    !! Estimates the frequency responses of a multiple-input, multiple-output
    !! (MIMO) system.
    real(real64), intent(in), dimension(:,:) :: x
        !! An N-by-P array containing the P inputs to the system.
    real(real64), intent(in), dimension(:,:) :: y
        !! An N-by-M array containing the M outputs from the system.
    real(real64), intent(in) :: fs
        !! The sampling frequency, in Hz.
    class(window), intent(in), optional, target :: win
        !! The window to apply to the data.  If nothing is supplied, no window
        !! is applied.
    integer(int32), intent(in), optional :: method
        !! Enter 1 to utilize an H1 estimator; else, enter 2 to utilize an
        !! H2 estimator.  The default is an H1 estimator.
        !!
        !! An H1 estimator is defined as the cross-spectrum of the input and
        !! response signals divided by the energy spectral density of the input.
        !! An H2 estimator is defined as the energy spectral density of the
        !! response divided by the cross-spectrum of the input and response
        !! signals.
        !!
        !! $$ H_{1} = \frac{P_{xy}}{P_{xx}} $$
        !!
        !! $$ H_{2} = \frac{P_{yy}}{P_{xy}} $$
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to provide 
        !! error handling. Possible errors and warning messages that may be 
        !! encountered are as follows.
        !!
        !! - DYN_MEMORY_ERROR: Occurs if there are issues allocating memory.
        !!
        !! - DYN_ARRAY_SIZE_ERROR: Occurs if x and y do not have the same number
        !!   of rows.
    type(mimo_frf) :: rst
        !! The resulting frequency response functions.

    ! Local Variables
    integer(int32) :: i, j, npts, m, p, nfreq, meth, flag
    real(real64) :: df
    class(window), pointer :: wptr
    type(rectangular_window), target :: defwin
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(x, 1)
    m = size(y, 2)
    p = size(x, 2)
    if (present(win)) then
        wptr => win
    else
        defwin%size = npts
        wptr => defwin
    end if
    if (present(method)) then
        if (method == 2) then
            meth = SPCTRM_H2_ESTIMATOR
        else
            meth = SPCTRM_H1_ESTIMATOR
        end if
    else
        meth = SPCTRM_H1_ESTIMATOR
    end if
    nfreq = compute_transform_length(wptr%size)
    allocate(rst%frequency(nfreq), stat = flag)
    if (flag == 0) allocate(rst%responses(nfreq, m, p))
    if (flag /= 0) then
        call report_memory_error("mimo_freqres", flag, errmgr)
        return
    end if

    ! Input Checking
    if (size(y, 1) /= npts) then
        call report_matrix_size_error("mimo_freqres", "y", npts, size(y, 2), &
            size(y, 1), size(y, 2), errmgr)
        return
    end if

    ! Compute the transfer functions for each possible combination
    do j = 1, p
        do i = 1, m
            rst%responses(:,i,j) = siso_transfer_function(wptr, &
                x(:,j), y(:,i), etype = meth, err = errmgr)
            if (errmgr%has_error_occurred()) return
        end do
    end do

    ! Compute the frequency vector
    df = frequency_bin_width(fs, wptr%size)
    rst%frequency = (/ (df * i, i = 0, nfreq - 1) /)
end function

! ******************************************************************************
! V1.0.6 ADDITIONS
! ------------------------------------------------------------------------------
! SEE: https://www.researchgate.net/publication/224619803_Reduction_of_structure-borne_noise_in_automobiles_by_multivariable_feedback
subroutine frf_accel_fit_fcn(xdata, mdl, rst, stop)
    !! The FRF fitting function for an accelerance FRF (acceleration-excited).
    real(real64), intent(in), dimension(:) :: xdata
        !! The independent variable data.
    real(real64), intent(in), dimension(:) :: mdl
        !! The model parameters.
    real(real64), intent(out), dimension(:) :: rst
        !! The model results.
    logical, intent(out) :: stop
        !! Stop the simulation?

    ! Local Variables
    integer(int32) :: i, n
    complex(real64) :: h

    ! Process:
    ! 
    ! The amplitude portion of the response is stored in the first "N" locations
    ! in the output with the phase portion (in radians) is stored in the
    ! second "N" locations.
    n = size(xdata) / 2
    do i = 1, n
        h = evaluate_accelerance_frf_model(mdl, xdata(i))
        rst(i) = abs(h)
        rst(i + n) = atan2(aimag(h), real(h))
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine frf_force_fit_fcn(xdata, mdl, rst, stop)
    !! The FRF fitting function for an force-excited FRF.
    real(real64), intent(in), dimension(:) :: xdata
        !! The independent variable data.
    real(real64), intent(in), dimension(:) :: mdl
        !! The model parameters.
    real(real64), intent(out), dimension(:) :: rst
        !! The model results.
    logical, intent(out) :: stop
        !! Stop the simulation?

    ! Local Variables
    integer(int32) :: i, n
    complex(real64) :: h

    ! Process:
    ! 
    ! The amplitude portion of the response is stored in the first "N" locations
    ! in the output with the phase portion (in radians) is stored in the
    ! second "N" locations.
    n = size(xdata) / 2
    do i = 1, n
        h = evaluate_receptance_frf_model(mdl, xdata(i))
        rst(i) = abs(h)
        rst(i + n) = atan2(aimag(h), real(h))
    end do
end subroutine

! ------------------------------------------------------------------------------
function fit_frf(mt, n, freq, rsp, maxp, minp, init, stats, alpha, controls, &
    settings, info, err) result(rst)
    use peaks
    !! Fits an experimentally obtained frequency response by model for either a
    !! receptance model:
    !!
    !! $$ H(\omega) = \sum_{i=1}^{n} \frac{A_{i}}{\omega_{ni}^{2} - 
    !! \omega^{2} + 2 j \zeta_{i} \omega_{ni} \omega} $$
    !!
    !! or an accelerance model:
    !!
    !! $$ H(\omega) = \sum_{i=1}^{n} \frac{-A_{i} \omega^{2}}{\omega_{ni}^{2} - 
    !! \omega^{2} + 2 j \zeta_{i} \omega_{ni} \omega} $$.
    !!
    !! Internally, the code uses a Levenberg-Marquardt solver to determine the
    !! parameters.  The initial guess for the solver is determined by a 
    !! peak finding algorithm used to locate the resonant modes in frequency.
    !! from this result, estimates for both the amplitude and natural frequency
    !! values are obtained.  The damping parameters are assumed to be equal
    !! for all modes and set to a default value of 0.1.
    integer(int32), intent(in) :: mt
        !! The excitation method.  The options are as follows.
        !!
        !! - FRF_ACCELERANCE_MODEL: Use an accelerance model.
        !!
        !! - FRF_RECEPTANCE_MODEL: Use a receptance model.
    integer(int32), intent(in) :: n
        !! The model order (# of resonant modes).
    real(real64), intent(in), dimension(:) :: freq
        !! An M-element array containing the excitation frequency values in 
        !! units of rad/s.
    complex(real64), intent(in), dimension(:) :: rsp
        !! An M-element array containing the frequency response to fit.
    real(real64), intent(in), dimension(:), optional :: maxp
        !! An optional 3*N-element array that can be used as upper limits on 
        !! the parameter values. If no upper limit is requested for a particular
        !! parameter, utilize a very large value. The internal default is to 
        !! utilize huge() as a value.
    real(real64), intent(in), dimension(:), optional :: minp
        !! An optional 3*N-element array that can be used as lower limits on 
        !! the parameter values. If no lower limit is requested for a particalar
        !! parameter, utilize a very large magnitude, but negative, value. The 
        !! internal default is to utilize -huge() as a value.
    real(real64), intent(in), dimension(:), optional :: init
        !! An optional 3*N-element array that, if supplied, provides an initial
        !! guess for each of the 3*N model parameters for the iterative solver.
        !! If supplied, this array replaces the peak finding algorithm for
        !! estimating an initial guess.
    type(regression_statistics), intent(out), dimension(:), optional :: stats
        !! An optional 3*N-element array that, if supplied, will be used to
        !! return statistics about the fit for each model parameter.
    real(real64), intent(in), optional :: alpha
        !! The significance level at which to evaluate the confidence intervals.
        !! The default value is 0.05 such that a 95% confidence interval is 
        !! calculated.
    type(iteration_controls), intent(in), optional :: controls
        !! An optional input providing custom iteration controls.
    type(lm_solver_options), intent(in), optional :: settings
        !! An optional input providing custom settings for the solver.
    type(convergence_info), intent(out), optional :: info
        !! An optional output that can be used to gain information about the 
        !! iterative solution and the nature of the convergence.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided 
        !! can be used to retrieve information relating to any errors 
        !! encountered during execution. If not provided, a default 
        !! implementation of the errors class is used internally to provide 
        !! error handling. Possible errors and warning messages that may be 
        !! encountered are as follows.
        !!
        !! - DYN_MEMORY_ERROR: Occurs if there are issues allocating memory.
        !!
        !! - DYN_ARRAY_SIZE_ERROR: Occurs if freq and rsp are not the same size.
        !!
        !! - DYN_UNDERDEFINED_PROBLEM_EROR: Occurs if the requested model 
        !!      order is too high for the number of data points available.
        !!
        !! - DYN_TOLERANCE_TOO_SMALL_ERROR: Occurs if the requested solver 
        !!      tolerance is too small to be practical for this problem.
        !!
        !! - DYN_TOO_FEW_ITERATIONS_ERROR: Occurs if convergence cannot be 
        !!      achieved in the allowed number of solver iterations.
    real(real64), allocatable, dimension(:) :: rst
        !! An array containing the model parameters stored as $$ \left[ A_{1}, 
        !! \omega_{n1}, \zeta_{1}, A_{2}, \omega_{n2}, \zeta_{2} ... \right] $$.

    ! Parameters
    real(real64), parameter :: zeta = 0.1d0

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    procedure(regression_function), pointer :: fcn
    integer(int32) :: i, npts, nparam, flag
    integer(int32), allocatable, dimension(:) :: maxinds, mininds
    real(real64) :: maxamp, minamp, amprange, delta
    real(real64), allocatable, dimension(:) :: x, y, maxvals, minvals, &
        ymod, resid
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    select case (mt)
    case (FRF_ACCELERANCE_MODEL)
        fcn => frf_accel_fit_fcn
    case (FRF_RECEPTANCE_MODEL)
        fcn => frf_force_fit_fcn
    case default
        call errmgr%report_error("fit_frf", "Invalid entry for the " // &
            "variable mt.  Must be either FRF_ACCELERANCE_MODEL or " // &
            "FRF_RECEPTANCE_MODEL.", DYN_INVALID_INPUT_ERROR)
        return
    end select
    npts = size(freq)
    nparam = 3 * n

    ! Input Checking
    if (size(rsp) /= npts) then
        call report_array_size_error("fit_frf", "rsp", npts, size(rsp), errmgr)
        return
    end if

    ! Memory Allocations
    allocate(rst(nparam), stat = flag)
    if (flag == 0) allocate( &
        x(2 * npts), &
        y(2 * npts), &
        ymod(2 * npts), &
        resid(2 * npts), &
        stat = flag)
    if (flag /= 0) then
        call report_memory_error("fit_frf", flag, errmgr)
        return
    end if

    ! Determine phase and amplitude terms, and store frequency values
    do i = 1, npts
        ! Store frequency values
        x(i) = freq(i)
        x(i + npts) = freq(i)

        ! Store amplitude and phase values
        y(i) = abs(rsp(i))
        y(i + npts) = atan2(aimag(rsp(i)), real(rsp(i)))

        ! Determine max and min amplitudes
        if (i == 1) then
            maxamp = y(i)
            minamp = y(i)
        else
            if (y(i) > maxamp) maxamp = y(i)
            if (y(i) < minamp) minamp = y(i)
        end if
    end do
    amprange = maxamp - minamp

    if (present(init)) then
        ! Check the array size
        if (size(init) /= nparam) then
            call report_array_size_error("fit_frf", "init", nparam, &
                size(init), errmgr)
            return
        end if

        ! Copy init to rst
        rst = init
    else
        ! Perform the peak location to determine an initial guess at parameters
        delta = 0.005d0 * amprange
        call peak_detect(y(1:npts), delta, maxinds, maxvals, mininds, minvals)
        do i = 1, min(n, size(maxvals))
            rst(3 * i - 2) = maxvals(i)         ! amplitude
            rst(3 * i - 1) = freq(maxinds(i))   ! frequency
            rst(3 * i) = zeta                   ! damping
        end do
        if (size(maxvals) < n) then
            ! The peak detection did not find enough peaks.
            if (size(maxvals) == 0) then
                ! No peaks found.  This is suspicious, but just use a random
                ! estimate to get started.  Maybe the solver will be able
                ! to sort it out.
                call random_number(rst)
            else
                ! Fill in the remaining parameters with the last set estimate
                do i = size(maxvals) + 1, n
                    rst(3 * i - 2) = rst(3 * (i - 1) - 2)
                    rst(3 * i - 1) = rst(3 * (i - 1) - 1)
                    rst(3 * i) = rst(3 * (i - 1))
                end do
            end if
        end if
    end if

    ! Fit the model
    call nonlinear_least_squares(fcn, x, y, rst, ymod, resid, maxp = maxp, &
        minp = minp, stats = stats, alpha = alpha, controls = controls, &
        settings = settings, info = info, err = errmgr)
end function

! ------------------------------------------------------------------------------
pure function evaluate_accelerance_frf_model_scalar(mdl, w) result(rst)
    !! Evaluates the specified accelerance FRF model.  The model is of
    !! the following form.
    !!
    !! $$ H(\omega) = \sum_{i=1}^{n} \frac{-A_{i} \omega^{2}}{\omega_{ni}^{2} - 
    !! \omega^{2} + 2 j \zeta_{i} \omega_{ni} \omega}  $$
    real(real64), intent(in), dimension(:) :: mdl
        !! The model parameter array.  The elements of the array are stored
        !! as $$ \left[ A_{1}, \omega_{n1}, \zeta_{1}, A_{2}, \omega_{n2}, 
        !! \zeta_{2} ... \right] $$.
    real(real64), intent(in) :: w
        !! The frequency value, in rad/s, at which to evaluate the model.
    complex(real64) :: rst
        !! The resulting frequency response function.

    ! Local Variables
    integer(int32) :: i, j, n

    ! Process
    j = 1
    n = size(mdl) / 3
    rst = (0.0d0, 0.0d0)
    do i = 1, n
        rst = rst + frf_accel_model_driver(mdl(j), mdl(j+1), mdl(j+2), w)
        j = j + 3
    end do
end function

! ----------
pure function evaluate_accelerance_frf_model_array(mdl, w) result(rst)
    !! Evaluates the specified accelerance FRF model.  The model is of
    !! the following form.
    !!
    !! $$ H(\omega) = \sum_{i=1}^{n} \frac{-A_{i} \omega^{2}}{\omega_{ni}^{2} - 
    !! \omega^{2} + 2 j \zeta_{i} \omega_{ni} \omega}  $$
    real(real64), intent(in), dimension(:) :: mdl
        !! The model parameter array.  The elements of the array are stored
        !! as $$ \left[ A_{1}, \omega_{n1}, \zeta_{1}, A_{2}, \omega_{n2}, 
        !! \zeta_{2} ... \right] $$.
    real(real64), intent(in), dimension(:) :: w
        !! The frequency value, in rad/s, at which to evaluate the model.
    complex(real64), allocatable, dimension(:) :: rst
        !! The resulting frequency response function.

    ! Local Variables
    integer(int32) :: i, n

    ! Process
    n = size(w)
    allocate(rst(n))
    do i = 1, n
        rst(i) = evaluate_accelerance_frf_model_scalar(mdl, w(i))
    end do
end function

! ----------
pure elemental function frf_accel_model_driver(A, wn, zeta, w) result(rst)
    !! Evaluates a single term of the accelerance FRF model.
    real(real64), intent(in) :: A
        !! The amplitude term.
    real(real64), intent(in) :: wn
        !! The natural frequency term.
    real(real64), intent(in) :: zeta
        !! The damping ratio term.
    real(real64), intent(in) :: w
        !! The excitation frequency.
    complex(real64) :: rst
        !! The result.

    ! Parameters
    complex(real64), parameter :: j = (0.0d0, 1.0d0)

    ! Process
    rst = -A * w**2 / (wn**2 - w**2 + 2.0d0 * j * zeta * wn * w)
end function

! ------------------------------------------------------------------------------
pure function evaluate_receptance_frf_model_scalar(mdl, w) result(rst)
    !! Evaluates the specified receptance FRF model.  The model is of
    !! the following form.
    !!
    !! $$ H(\omega) = \sum_{i=1}^{n} \frac{A_{i}}{\omega_{ni}^{2} - 
    !! \omega^{2} + 2 j \zeta_{i} \omega_{ni} \omega}  $$
    real(real64), intent(in), dimension(:) :: mdl
        !! The model parameter array.  The elements of the array are stored
        !! as $$ \left[ A_{1}, \omega_{n1}, \zeta_{1}, A_{2}, \omega_{n2}, 
        !! \zeta_{2} ... \right] $$.
    real(real64), intent(in) :: w
        !! The frequency value, in rad/s, at which to evaluate the model.
    complex(real64) :: rst
        !! The resulting frequency response function.

    ! Local Variables
    integer(int32) :: i, j, n

    ! Process
    j = 1
    n = size(mdl) / 3
    rst = (0.0d0, 0.0d0)
    do i = 1, n
        rst = rst + frf_receptance_model_driver(mdl(j), mdl(j+1), mdl(j+2), w)
        j = j + 3
    end do
end function

! ----------
pure function evaluate_receptance_frf_model_array(mdl, w) result(rst)
    !! Evaluates the specified receptance FRF model.  The model is of
    !! the following form.
    !!
    !! $$ H(\omega) = \sum_{i=1}^{n} \frac{A_{i}}{\omega_{ni}^{2} - 
    !! \omega^{2} + 2 j \zeta_{i} \omega_{ni} \omega}  $$
    real(real64), intent(in), dimension(:) :: mdl
        !! The model parameter array.  The elements of the array are stored
        !! as $$ \left[ A_{1}, \omega_{n1}, \zeta_{1}, A_{2}, \omega_{n2}, 
        !! \zeta_{2} ... \right] $$.
    real(real64), intent(in), dimension(:) :: w
        !! The frequency value, in rad/s, at which to evaluate the model.
    complex(real64), allocatable, dimension(:) :: rst
        !! The resulting frequency response function.

    ! Local Variables
    integer(int32) :: i, n

    ! Process
    n = size(w)
    allocate(rst(n))
    do i = 1, n
        rst(i) = evaluate_receptance_frf_model_scalar(mdl, w(i))
    end do
end function

! ----------
pure elemental function frf_receptance_model_driver(A, wn, zeta, w) result(rst)
    !! Evaluates a single term of the receptance FRF model.
    real(real64), intent(in) :: A
        !! The amplitude term.
    real(real64), intent(in) :: wn
        !! The natural frequency term.
    real(real64), intent(in) :: zeta
        !! The damping ratio term.
    real(real64), intent(in) :: w
        !! The excitation frequency.
    complex(real64) :: rst
        !! The result.

    ! Parameters
    complex(real64), parameter :: j = (0.0d0, 1.0d0)

    ! Process
    rst = A / (wn**2 - w**2 + 2.0d0 * j * zeta * wn * w)
end function

! ------------------------------------------------------------------------------
end module
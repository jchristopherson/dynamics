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
module dynamics_frequency_response
    use iso_fortran_env
    use dynamics_error_handling
    use dynamics_modal_analysis
    use spectrum
    use fstats
    use dynamics_helper
    use lapack, only : ZGELSY
    implicit none
    private
    public :: modal_excite
    public :: frf
    public :: mimo_frf
    public :: frequency_response
    public :: evaluate_accelerance_frf_model
    public :: evaluate_receptance_frf_model
    public :: fit_frf
    public :: FRF_ACCELERANCE_MODEL
    public :: FRF_RECEPTANCE_MODEL
    public :: regression_statistics
    public :: iteration_controls
    public :: lm_solver_options
    public :: convergence_info
    public :: dynamic_stiffness

    interface
        subroutine modal_excite(freq, frc, args)
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
            class(*), intent(inout), optional :: args
                !! An optional argument that can be used to communicate with
                !! the outside world.
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
            !! A P-element array containing the frequency values at which the 
            !! FRF is provided.  The units of this array are the same as the
            !! units of the frequency values passed to the routine used to 
            !! compute the frequency response.
        complex(real64), allocatable, dimension(:,:,:) :: responses
            !! An N-by-M-by-P array containing the N frequency response 
            !! functions for each of the M inputs corresponding to each of 
            !! the P frequency points.
    end type

    interface frequency_response
        !! Computes the frequency response functions for a system of ODE's.
        module procedure :: frf_modal_prop_damp
        module procedure :: frf_modal_prop_damp_sparse
        module procedure :: frf_modal_prop_damp_2
        module procedure :: frf_modal_prop_damp_sparse_2
        module procedure :: frf_general_damp_1
        module procedure :: frf_general_damp_2
        module procedure :: siso_freqres
        module procedure :: mimo_freqres
    end interface

    interface evaluate_accelerance_frf_model
        module procedure :: evaluate_accelerance_frf_model_scalar
        module procedure :: evaluate_accelerance_frf_model_array
    end interface

    interface evaluate_receptance_frf_model
        module procedure :: evaluate_receptance_frf_model_scalar
        module procedure :: evaluate_receptance_frf_model_array
    end interface

    interface dynamic_stiffness
        module procedure :: dynamic_stiffness_dense
    end interface

! ------------------------------------------------------------------------------
    integer(int32), parameter :: FRF_ACCELERANCE_MODEL = 1
        !! Defines an accelerance frequency response model.
    integer(int32), parameter :: FRF_RECEPTANCE_MODEL = 2
        !! Defines a receptance frequency response model.

contains
! ------------------------------------------------------------------------------
    function frf_modal_prop_damp(mass, stiff, alpha, beta, freq, frc, &
        modes, modeshapes, args) result(rst)
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
        class(*), intent(inout), optional :: args
            !! An optional argument that can be used to communicate with
            !! the outside world.
        type(frf) :: rst
            !! The resulting frequency responses.

        ! Parameters
        complex(real64), parameter :: j = (0.0d0, 1.0d0)
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        ! Local Variables
        integer(int32) :: i, m, n
        complex(real64) :: s
        real(real64), allocatable, dimension(:) :: lambda, zeta
        complex(real64), allocatable, dimension(:) :: vals, q, f, u
        complex(real64), allocatable, dimension(:,:) :: vecs
        
        ! Initialization
        m = size(freq)
        n = size(mass, 1)

        ! Input Checking
        if (n < 1) error stop DYN_INVALID_INPUT_ERROR
        if (size(mass, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR
        if (size(stiff,1) /= size(stiff, 2)) error stop DYN_MATRIX_SIZE_ERROR
        if (size(stiff, 1) /= n .or. size(stiff, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR
        if (.not.is_symmetric(mass) .or. .not.is_symmetric(stiff)) &
            error stop DYN_INVALID_INPUT_ERROR
        if (.not.(alpha >= 0.0d0) .or. .not.(beta >= 0.0d0)) &
            error stop DYN_INVALID_INPUT_ERROR
        if (.not.associated(frc)) error stop DYN_NULL_POINTER_ERROR

        ! Memory allocations
        allocate(zeta(n))
        allocate(q(n))
        allocate(vals(n))
        allocate(f(n), source = zero)
        allocate(u(n))
        allocate(vecs(n, n))
        allocate(rst%responses(m, n))
        allocate(rst%frequency(m), source = freq)

        ! Compute the eigenvalues and eigenvectors
        call eigen(stiff, mass, vals, rvecs = vecs)
        allocate(lambda(n), source = real(vals))
        if (any(lambda <= 0.0d0)) error stop DYN_INVALID_INPUT_ERROR

        ! Compute the damping terms
        zeta = compute_modal_damping(lambda, alpha, beta)

        ! Compute each transfer function
        do i = 1, m
            call frc(freq(i), f, args)
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
            allocate(modes(n), source = sqrt(real(vals)))
        end if

        if (present(modeshapes)) then
            allocate(modeshapes(n, n), source = real(vecs))
        end if
    end function

! ------------------------------------------------------------------------------
    function frf_modal_prop_damp_sparse(mass, stiff, alpha, beta, nmodes, &
        freq, frc, modes, modeshapes, args) result(rst)
        !! Computes a modal-truncated frequency response for a system with
        !! proportional damping using CSR sparse mass and stiffness matrices.
        !! The damping matrix is defined by \(C=\alpha M+\beta K\).
        use dynamics_error_handling
        use linalg, only : csr_matrix, matmul, size
        type(csr_matrix), intent(in) :: mass
            !! The N-by-N symmetric positive-definite mass matrix.
        type(csr_matrix), intent(in) :: stiff
            !! The N-by-N symmetric stiffness matrix.
        real(real64), intent(in) :: alpha
            !! The mass damping factor, \(\alpha\).
        real(real64), intent(in) :: beta
            !! The stiffness damping factor, \(\beta\).
        integer(int32), intent(in) :: nmodes
            !! The number of lowest-frequency modes to retain.  This value
            !! must be greater than zero and less than N.
        real(real64), intent(in), dimension(:) :: freq
            !! An M-element array of frequency values in units of rad/s.
        procedure(modal_excite), pointer, intent(in) :: frc
            !! A pointer to the physical forcing function.
        real(real64), intent(out), allocatable, optional, dimension(:) :: modes
            !! An optional NMODES-element array containing the retained modal
            !! frequencies in units of rad/s.
        real(real64), intent(out), allocatable, optional, dimension(:,:) :: &
            modeshapes
            !! An optional N-by-NMODES matrix containing the mass-normalized
            !! retained mode shapes.
        class(*), intent(inout), optional :: args
            !! An optional argument passed to the forcing function.
        type(frf) :: rst
            !! The modal-truncated frequency responses.

        complex(real64), parameter :: j = (0.0d0, 1.0d0)
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)

        integer(int32) :: i, imode, m, n
        real(real64) :: modal_mass
        complex(real64) :: s
        real(real64), allocatable, dimension(:) :: mass_vec, modal_freqs, zeta
        real(real64), allocatable, dimension(:,:) :: vecs
        complex(real64), allocatable, dimension(:) :: f, q, u

        m = size(freq)
        n = size(mass, 1)

        if (.not.(alpha >= 0.0d0) .or. .not.(beta >= 0.0d0)) &
            error stop DYN_INVALID_INPUT_ERROR
        if (.not.associated(frc)) error stop DYN_NULL_POINTER_ERROR

        call modal_response(mass, stiff, nmodes, modal_freqs, vecs)

        allocate(mass_vec(n))
        do imode = 1, nmodes
            mass_vec = matmul(mass, vecs(:,imode))
            modal_mass = dot_product(vecs(:,imode), mass_vec)
            if (.not.(modal_mass > 0.0d0)) &
                error stop DYN_INVALID_INPUT_ERROR
            vecs(:,imode) = vecs(:,imode) / sqrt(modal_mass)
        end do

        allocate(zeta(nmodes), source = compute_modal_damping( &
            modal_freqs**2, alpha, beta))
        allocate(f(n), source = zero)
        allocate(q(nmodes), source = zero)
        allocate(u(nmodes), source = zero)
        allocate(rst%responses(m, n), source = zero)
        allocate(rst%frequency(m), source = freq)

        do i = 1, m
            call frc(freq(i), f, args)
            do imode = 1, nmodes
                u(imode) = sum(vecs(:,imode) * f)
            end do
            s = j * freq(i)
            q = u / (s**2 + 2.0d0 * zeta * modal_freqs * s + &
                modal_freqs**2)
            do imode = 1, nmodes
                rst%responses(i,:) = rst%responses(i,:) + &
                    vecs(:,imode) * q(imode)
            end do
        end do

        if (present(modes)) then
            allocate(modes(nmodes), source = modal_freqs)
        end if
        if (present(modeshapes)) then
            allocate(modeshapes(n, nmodes), source = vecs)
        end if
    end function

! ------------------------------------------------------------------------------
    function frf_modal_prop_damp_2(mass, stiff, alpha, beta, nfreq, freq1, &
        freq2, frc, modes, modeshapes, args) result(rst)
        !! Computes the frequency response functions for a 
        !! multi-degree-of-freedom system that uses proportional damping such
        !! that the damping matrix \( C \) is related to the stiffness an mass
        !! matrices by proportional damping coefficients \( \alpha \) and
        !! In modal coordinates, each mode has denominator
        !! $$ s^2+2\zeta_i\omega_i s+\omega_i^2, $$
        !! and the physical response is reconstructed from the mode shapes.
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
        class(*), intent(inout), optional :: args
            !! An optional argument that can be used to communicate with
            !! the outside world.
        type(frf) :: rst
            !! The resulting frequency responses.

        ! Local Variables
        integer(int32) :: i, flag
        real(real64) :: df
        real(real64), allocatable, dimension(:) :: freq

        ! Input Checking
        if (abs(freq1 - freq2) < sqrt(epsilon(freq1))) error stop DYN_INVALID_INPUT_ERROR
        if (nfreq < 2) error stop DYN_INVALID_INPUT_ERROR

        ! Process
        df = (freq2 - freq1) / (nfreq - 1.0d0)
        allocate(freq(nfreq))
        freq = (/ (df * i + freq1, i = 0, nfreq - 1) /)
        rst = frequency_response(mass, stiff, alpha, beta, freq, frc, modes, &
            modeshapes, args = args)
    end function

! ------------------------------------------------------------------------------
    function frf_modal_prop_damp_sparse_2(mass, stiff, alpha, beta, nmodes, &
        nfreq, freq1, freq2, frc, modes, modeshapes, args) result(rst)
        !! Computes a modal-truncated frequency response for a system with
        !! proportional damping using CSR sparse mass and stiffness matrices.
        !! The damping matrix is defined by \(C=\alpha M+\beta K\).
        use dynamics_error_handling
        use linalg, only : csr_matrix, matmul, size
        type(csr_matrix), intent(in) :: mass
            !! The N-by-N symmetric positive-definite mass matrix.
        type(csr_matrix), intent(in) :: stiff
            !! The N-by-N symmetric stiffness matrix.
        real(real64), intent(in) :: alpha
            !! The mass damping factor, \(\alpha\).
        real(real64), intent(in) :: beta
            !! The stiffness damping factor, \(\beta\).
        integer(int32), intent(in) :: nmodes
            !! The number of lowest-frequency modes to retain.  This value
            !! must be greater than zero and less than N.
        integer(int32), intent(in) :: nfreq
            !! The number of frequency values to analyze.  This value must be
            !! at least 2.
        real(real64), intent(in) :: freq1
            !! The starting frequency, in units of rad/s.
        real(real64), intent(in) :: freq2
            !! The ending frequency, in units of rad/s.
        procedure(modal_excite), pointer, intent(in) :: frc
            !! A pointer to the physical forcing function.
        real(real64), intent(out), allocatable, optional, dimension(:) :: modes
            !! An optional NMODES-element array containing the retained modal
            !! frequencies in units of rad/s.
        real(real64), intent(out), allocatable, optional, dimension(:,:) :: &
            modeshapes
            !! An optional N-by-NMODES matrix containing the mass-normalized
            !! retained mode shapes.
        class(*), intent(inout), optional :: args
            !! An optional argument passed to the forcing function.
        type(frf) :: rst
            !! The modal-truncated frequency responses.

        ! Local Variables
        integer(int32) :: i
        real(real64) :: df
        real(real64), allocatable, dimension(:) :: freq

        ! Input Checking
        if (abs(freq1 - freq2) < sqrt(epsilon(freq1))) error stop DYN_INVALID_INPUT_ERROR
        if (nfreq < 2) error stop DYN_INVALID_INPUT_ERROR

        ! Process
        df = (freq2 - freq1) / (nfreq - 1.0d0)
        allocate(freq(nfreq))
        freq = (/ (df * i + freq1, i = 0, nfreq - 1) /)
        rst = frequency_response(mass, stiff, alpha, beta, nmodes, freq, frc, &
            modes, modeshapes, args = args)
    end function

! ******************************************************************************
! VERSION 1.0.5 ADDITIONS
! ------------------------------------------------------------------------------
function siso_freqres(x, y, fs, win, method) result(rst)
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
    type(frf) :: rst
        !! The resulting frequency response function.

    ! Local Variables
    integer(int32) :: i, npts, nfreq, meth
    real(real64) :: df
    class(window), pointer :: wptr
    type(rectangular_window), target :: defwin
    
    ! Input Checking
    npts = size(x)
    if (npts < 2 .or. .not.(fs > 0.0d0)) error stop DYN_INVALID_INPUT_ERROR
    if (size(y) /= npts) error stop DYN_ARRAY_SIZE_ERROR
    if (present(win)) then
        wptr => win
    else
        defwin%size = npts
        wptr => defwin
    end if
    if (present(method)) then
        if (method /= 1 .and. method /= 2) error stop DYN_INVALID_INPUT_ERROR
        if (method == 2) then
            meth = SPCTRM_H2_ESTIMATOR
        else
            meth = SPCTRM_H1_ESTIMATOR
        end if
    else
        meth = SPCTRM_H1_ESTIMATOR
    end if
    if (wptr%size < 2) error stop DYN_ARRAY_SIZE_ERROR
    nfreq = compute_transform_length(wptr%size)
    allocate(rst%frequency(nfreq))
    allocate(rst%responses(nfreq, 1))

    ! Compute the transfer function
    rst%responses(:,1) = siso_transfer_function(wptr, x, y, etype = meth)

    ! Compute the frequency vector
    df = frequency_bin_width(fs, wptr%size)
    rst%frequency = (/ (df * i, i = 0, nfreq - 1) /)
end function

! ------------------------------------------------------------------------------
function mimo_freqres(x, y, fs, win, method) result(rst)
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
    type(mimo_frf) :: rst
        !! The resulting frequency response functions.

    ! Local Variables
    integer(int32) :: i, j, npts, m, p, nfreq, meth
    real(real64) :: df
    class(window), pointer :: wptr
    type(rectangular_window), target :: defwin
    
    ! Input Checking
    npts = size(x, 1)
    m = size(y, 2)
    p = size(x, 2)
    if (npts < 2 .or. m < 1 .or. p < 1 .or. .not.(fs > 0.0d0)) &
        error stop DYN_INVALID_INPUT_ERROR
    if (size(y, 1) /= npts) error stop DYN_MATRIX_SIZE_ERROR
    if (present(win)) then
        wptr => win
    else
        defwin%size = npts
        wptr => defwin
    end if
    if (present(method)) then
        if (method /= 1 .and. method /= 2) error stop DYN_INVALID_INPUT_ERROR
        if (method == 2) then
            meth = SPCTRM_H2_ESTIMATOR
        else
            meth = SPCTRM_H1_ESTIMATOR
        end if
    else
        meth = SPCTRM_H1_ESTIMATOR
    end if
    if (wptr%size < 2) error stop DYN_ARRAY_SIZE_ERROR
    nfreq = compute_transform_length(wptr%size)
    allocate(rst%frequency(nfreq))

    ! Compute the transfer functions for each possible combination
    rst%responses = mimo_transfer_function(wptr, x, y, meth)

    ! Compute the frequency vector
    df = frequency_bin_width(fs, wptr%size)
    rst%frequency = (/ (df * i, i = 0, nfreq - 1) /)
end function

! ******************************************************************************
! V1.0.6 ADDITIONS
! ------------------------------------------------------------------------------
! SEE: https://www.researchgate.net/publication/224619803_Reduction_of_structure-borne_noise_in_automobiles_by_multivariable_feedback
subroutine frf_accel_fit_fcn(xdata, mdl, rst, stop, args)
    !! The FRF fitting function for an accelerance FRF (acceleration-excited).
    real(real64), intent(in), dimension(:) :: xdata
        !! The independent variable data.
    real(real64), intent(in), dimension(:) :: mdl
        !! The model parameters.
    real(real64), intent(out), dimension(:) :: rst
        !! The model results.
    logical, intent(out) :: stop
        !! Stop the simulation?
    class(*), intent(inout), optional :: args
        !! Optional arguments from the calling code.

    ! Local Variables
    integer(int32) :: i, n
    complex(real64) :: h

    ! Process:
    ! 
    ! The amplitude portion of the response is stored in the first "N" locations
    ! in the output with the phase portion (in radians) is stored in the
    ! second "N" locations.
    stop = .false.
    n = size(xdata) / 2
    do i = 1, n
        h = evaluate_accelerance_frf_model(mdl, xdata(i))
        rst(i) = abs(h)
        rst(i + n) = atan2(aimag(h), real(h))
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine frf_force_fit_fcn(xdata, mdl, rst, stop, args)
    !! The FRF fitting function for an force-excited FRF.
    real(real64), intent(in), dimension(:) :: xdata
        !! The independent variable data.
    real(real64), intent(in), dimension(:) :: mdl
        !! The model parameters.
    real(real64), intent(out), dimension(:) :: rst
        !! The model results.
    logical, intent(out) :: stop
        !! Stop the simulation?
    class(*), intent(inout), optional :: args
        !! Optional arguments from the calling code.

    ! Local Variables
    integer(int32) :: i, n
    complex(real64) :: h

    ! Process:
    ! 
    ! The amplitude portion of the response is stored in the first "N" locations
    ! in the output with the phase portion (in radians) is stored in the
    ! second "N" locations.
    stop = .false.
    n = size(xdata) / 2
    do i = 1, n
        h = evaluate_receptance_frf_model(mdl, xdata(i))
        rst(i) = abs(h)
        rst(i + n) = atan2(aimag(h), real(h))
    end do
end subroutine

! ------------------------------------------------------------------------------
function fit_frf(mt, n, freq, rsp, maxp, minp, init, stats, alpha, controls, &
    settings, info) result(rst)
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
    real(real64), allocatable, dimension(:) :: rst
        !! An array containing the model parameters stored as $$ \left[ A_{1}, 
        !! \omega_{n1}, \zeta_{1}, A_{2}, \omega_{n2}, \zeta_{2} ... \right] $$.

    ! Parameters
    real(real64), parameter :: zeta = 0.1d0

    ! Local Variables
    procedure(regression_function), pointer :: fcn
    integer(int32) :: i, npts, nparam
    integer(int32), allocatable, dimension(:) :: maxinds, mininds
    real(real64) :: maxamp, minamp, amprange, delta
    real(real64), allocatable, dimension(:) :: x, y, maxvals, minvals, &
        ymod, resid
    
    ! Initialization
    select case (mt)
    case (FRF_ACCELERANCE_MODEL)
        fcn => frf_accel_fit_fcn
    case (FRF_RECEPTANCE_MODEL)
        fcn => frf_force_fit_fcn
    case default
        error stop DYN_INVALID_INPUT_ERROR
    end select
    npts = size(freq)
    nparam = 3 * n

    ! Input Checking
    if (n < 1 .or. npts < 2) error stop DYN_INVALID_INPUT_ERROR
    if (size(rsp) /= npts) error stop DYN_ARRAY_SIZE_ERROR
    if (present(maxp)) then
        if (size(maxp) /= nparam) error stop DYN_ARRAY_SIZE_ERROR
    end if
    if (present(minp)) then
        if (size(minp) /= nparam) error stop DYN_ARRAY_SIZE_ERROR
    end if
    if (present(init)) then
        if (size(init) /= nparam) error stop DYN_ARRAY_SIZE_ERROR
    end if
    if (present(stats)) then
        if (size(stats) /= nparam) error stop DYN_ARRAY_SIZE_ERROR
    end if
    if (present(maxp) .and. present(minp)) then
        if (any(minp > maxp)) error stop DYN_INVALID_INPUT_ERROR
    end if

    ! Memory Allocations
    allocate( &
        rst(nparam), &
        x(2 * npts), &
        y(2 * npts), &
        ymod(2 * npts), &
        resid(2 * npts) &
    )

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
                ! No peaks found.  This is suspicious, but use a deterministic
                ! estimate to ensure predictable behavior.
                do i = 1, n
                    rst(3 * i - 2) = maxamp
                    rst(3 * i - 1) = freq(max(1, min(npts, (i * npts) / (n + 1))))
                    rst(3 * i) = zeta
                end do
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
        settings = settings, info = info)
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

! ******************************************************************************
! V1.9 ADDITIONS
! ------------------------------------------------------------------------------
    pure subroutine dynamic_stiffness_dense(omega, mass, damp, stiff, dyn_stiff)
        !! Computes the dynamic stiffness matrix at the specified frequency 
        !! such that /( K_{dyn}(\omega) = K - \omega^{2} M + j \omega C /).
        real(real64), intent(in) :: omega
            !! The frequency, in rad/s.
        real(real64), intent(in), dimension(:,:) :: mass
            !! The N-by-N mass matrix.
        real(real64), intent(in), dimension(:,:) :: damp
            !! The N-by-N damping matrix.
        real(real64), intent(in), dimension(:,:) :: stiff
            !! The N-by-N stiffness matrix.
        complex(real64), intent(out), dimension(:,:) :: dyn_stiff
            !! The N-by-N dynamic stiffness matrix.

        ! Parameters
        complex(real64), parameter :: j = (0.0d0, 1.0d0)

        ! Local Variables
        integer(int32) :: n

        ! Input Checking
        n = size(mass, 1)
        if (size(mass, 2) /= n) error stop DYN_NONSQUARE_MATRIX_ERROR
        if (size(damp, 1) /= n .or. size(damp, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR
        if (size(stiff, 1) /= n .or. size(stiff, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR
        if (size(dyn_stiff, 1) /= n .or. size(dyn_stiff, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR

        ! Process
        dyn_stiff = stiff - omega**2 * mass + j * omega * damp
    end subroutine

! ------------------------------------------------------------------------------
    function frf_general_damp_1(mass, damp, stiff, freq, frc, ranks, args) result(rst)
        !! Computes the frequency response functions for a multi-degree-of-freedom
        !! system that has a general damping matrix, and is not necessarily 
        !! symmetric.  The problem is treated as the solution to the linear
        !! system \( \left( K - \omega^{2} M + j \omega C \right) H(\omega) = 
        !! F(\omega) \).
        real(real64), intent(in), dimension(:,:) :: mass
            !! The N-by-N mass matrix.
        real(real64), intent(in), dimension(:,:) :: damp
            !! The N-by-N damping matrix.
        real(real64), intent(in), dimension(:,:) :: stiff
            !! The N-by-N stiffness matrix.
        real(real64), intent(in), dimension(:) :: freq
            !! An M-element array of frequency values at which to evaluate the
            !! frequency response functions, in units of rad/s.
        procedure(modal_excite), pointer, intent(in) :: frc
            !! A pointer to a routine used to compute the modal forcing 
            !! function.
        integer(int32), intent(out), optional, dimension(:) :: ranks
            !! Provides information on the rank of the dynamic stiffness matrix
            !! for each frequency.  If provided, this array must be the same
            !! length as freq.
        class(*), intent(inout), optional :: args
            !! An optional argument that can be used to communicate with
            !! the outside world.
        type(frf) :: rst
            !! The resulting frequency responses.

        ! Parallel Threshold
        integer(int32), parameter :: parallel_threshold = 1000

        ! Local Variables
        integer(int32) :: m, n, check

        ! Initialization
        m = size(freq)
        n = size(mass, 1)
        check = m * n**3

        ! Determine whether this problem should be solved in parallel
        if (check > parallel_threshold) then
            rst = frf_general_damp_parallel(mass, damp, stiff, freq, frc, ranks, args)
        else
            rst = frf_general_damp_serial(mass, damp, stiff, freq, frc, ranks, args)
        end if
    end function

! --------------------
    function frf_general_damp_serial(mass, damp, stiff, freq, frc, ranks, args) result(rst)
        !! Computes the frequency response functions for a multi-degree-of-freedom
        !! system that has a general damping matrix, and is not necessarily 
        !! symmetric.  The problem is treated as the solution to the linear
        !! system \( \left( K - \omega^{2} M + j \omega C \right) H(\omega) = 
        !! F(\omega) \).
        real(real64), intent(in), dimension(:,:) :: mass
            !! The N-by-N mass matrix.
        real(real64), intent(in), dimension(:,:) :: damp
            !! The N-by-N damping matrix.
        real(real64), intent(in), dimension(:,:) :: stiff
            !! The N-by-N stiffness matrix.
        real(real64), intent(in), dimension(:) :: freq
            !! An M-element array of frequency values at which to evaluate the
            !! frequency response functions, in units of rad/s.
        procedure(modal_excite), pointer, intent(in) :: frc
            !! A pointer to a routine used to compute the modal forcing 
            !! function.
        integer(int32), intent(out), optional, dimension(:) :: ranks
            !! Provides information on the rank of the dynamic stiffness matrix
            !! for each frequency.  If provided, this array must be the same
            !! length as freq.
        class(*), intent(inout), optional :: args
            !! An optional argument that can be used to communicate with
            !! the outside world.
        type(frf) :: rst
            !! The resulting frequency responses.

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: j = (0.0d0, 1.0d0)

        ! Local Variables
        integer(int32) :: i, m, n, lwork, lrwork, rnk, info
        integer(int32), allocatable, dimension(:) :: jpvt
        real(real64) :: rcond
        real(real64), allocatable, dimension(:) :: rwork
        complex(real64), allocatable, dimension(:) :: work
        complex(real64), allocatable, dimension(:,:) :: K_dyn
        complex(real64) :: dummy(1), temp(1)

        ! Input Checking
        m = size(freq)
        n = size(mass, 1)
        if (size(mass, 2) /= n) error stop DYN_NONSQUARE_MATRIX_ERROR
        if (size(damp, 1) /= n .or. size(damp, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR
        if (size(stiff, 1) /= n .or. size(stiff, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR
        if (present(ranks)) then
            if (size(ranks) /= m) error stop DYN_ARRAY_SIZE_ERROR
        end if

        ! Initialization
        lrwork = 2 * n
        rcond = epsilon(rcond)

        ! Memory Allocations
        allocate(rst%frequency(m), source = freq)
        allocate( &
            rst%responses(m, n), &
            jpvt(n), &
            rwork(lrwork), &
            K_dyn(n, n) &
        )

        ! Determine an appropriate workspace
        call ZGELSY(n, n, 1, K_dyn, n, dummy, n, jpvt, rcond, rnk, temp, &
            -1, rwork, info)
        lwork = int(temp(1), kind = int32)
        allocate(work(lwork))

        ! Loop over each frequency and solve the linear system
        do i = 1, m
            ! Evaluate the forcing function
            call frc(freq(i), rst%responses(i,:), args)

            ! Evaluate the dynamic stiffness
            call dynamic_stiffness(freq(i), mass, damp, stiff, K_dyn)

            ! Solve the linear system
            jpvt = 0
            call ZGELSY(n, n, 1, K_dyn, n, rst%responses(i,:), n, jpvt, rcond, &
                rnk, work, lwork, rwork, info)

            ! Store the rank?
            if (present(ranks)) then
                ranks(i) = rnk
            end if
        end do
    end function

! --------------------
    function frf_general_damp_parallel(mass, damp, stiff, freq, frc, ranks, &
        args) result(rst)
        !! Computes the frequency response functions for a multi-degree-of-freedom
        !! system that has a general damping matrix, and is not necessarily 
        !! symmetric.  The problem is treated as the solution to the linear
        !! system \( \left( K - \omega^{2} M + j \omega C \right) H(\omega) = 
        !! F(\omega) \).
        real(real64), intent(in), dimension(:,:) :: mass
            !! The N-by-N mass matrix.
        real(real64), intent(in), dimension(:,:) :: damp
            !! The N-by-N damping matrix.
        real(real64), intent(in), dimension(:,:) :: stiff
            !! The N-by-N stiffness matrix.
        real(real64), intent(in), dimension(:) :: freq
            !! An M-element array of frequency values at which to evaluate the
            !! frequency response functions, in units of rad/s.
        procedure(modal_excite), pointer, intent(in) :: frc
            !! A pointer to a routine used to compute the modal forcing 
            !! function.
        integer(int32), intent(out), optional, dimension(:) :: ranks
            !! Provides information on the rank of the dynamic stiffness matrix
            !! for each frequency.  If provided, this array must be the same
            !! length as freq.
        class(*), intent(inout), optional :: args
            !! An optional argument that can be used to communicate with
            !! the outside world.
        type(frf) :: rst
            !! The resulting frequency responses.

        ! Parameters
        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: j = (0.0d0, 1.0d0)

        ! Local Variables
        logical :: return_rank
        integer(int32) :: i, m, n, lwork, lrwork, rnk, info, idummy(1)
        integer(int32), allocatable, dimension(:) :: jpvt
        real(real64) :: rcond, rdummy(1)
        real(real64), allocatable, dimension(:) :: rwork
        complex(real64), allocatable, dimension(:) :: work, f
        complex(real64), allocatable, dimension(:,:) :: K_dyn
        complex(real64) :: dummy(1), temp(1)

        ! Input Checking
        m = size(freq)
        n = size(mass, 1)
        if (size(mass, 2) /= n) error stop DYN_NONSQUARE_MATRIX_ERROR
        if (size(damp, 1) /= n .or. size(damp, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR
        if (size(stiff, 1) /= n .or. size(stiff, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR
        if (present(ranks)) then
            if (size(ranks) /= m) error stop DYN_ARRAY_SIZE_ERROR
        end if

        ! Initialization
        lrwork = 2 * n
        rcond = epsilon(rcond)
        return_rank = present(ranks)

        ! Memory Allocations
        allocate(rst%frequency(m), source = freq)
        allocate(rst%responses(m, n))

        ! Determine an appropriate workspace
        call ZGELSY(n, n, 1, dummy, n, dummy, n, idummy, rcond, rnk, temp, &
            -1, rdummy, info)
        lwork = int(temp(1), kind = int32)
        allocate(work(lwork))

        ! Process
        if (present(args)) then
            !$omp parallel do private(f, jpvt, K_dyn, work, rwork, rnk, info, args)
            do i = 1, m
                ! Memory Allocations
                if (.not.allocated(f)) allocate(f(n))
                if (.not.allocated(jpvt)) allocate(jpvt(n))
                if (.not.allocated(K_dyn)) allocate(K_dyn(n, n))
                if (.not.allocated(work)) allocate(work(lwork))
                if (.not.allocated(rwork)) allocate(rwork(lrwork))

                ! Evaluate the forcing function
                call frc(freq(i), f, args)

                ! Evaluate the dynamic stiffness
                call dynamic_stiffness(freq(i), mass, damp, stiff, K_dyn)

                ! Solve the linear system
                jpvt = 0
                call ZGELSY(n, n, 1, K_dyn, n, f, n, jpvt, rcond, rnk, work, &
                    lwork, rwork, info)
                
                ! Store the output
                rst%responses(i,:) = f
                if (return_rank) then
                    ranks(i) = rnk
                end if
            end do
            !$omp end parallel do
        else
            !$omp parallel do private(f, jpvt, K_dyn, work, rwork, rnk, info)
            do i = 1, m
                ! Memory Allocations
                if (.not.allocated(f)) allocate(f(n))
                if (.not.allocated(jpvt)) allocate(jpvt(n))
                if (.not.allocated(K_dyn)) allocate(K_dyn(n, n))
                if (.not.allocated(work)) allocate(work(lwork))
                if (.not.allocated(rwork)) allocate(rwork(lrwork))

                ! Evaluate the forcing function
                call frc(freq(i), f)

                ! Evaluate the dynamic stiffness
                call dynamic_stiffness(freq(i), mass, damp, stiff, K_dyn)

                ! Solve the linear system
                jpvt = 0
                call ZGELSY(n, n, 1, K_dyn, n, f, n, jpvt, rcond, rnk, work, &
                    lwork, rwork, info)
                
                ! Store the output
                rst%responses(i,:) = f
                if (return_rank) then
                    ranks(i) = rnk
                end if
            end do
            !$omp end parallel do
        end if
    end function

! ------------------------------------------------------------------------------
    function frf_general_damp_2(mass, damp, stiff, nfreq, freq1, freq2, frc, &
        ranks, args) result(rst)
        !! Computes the frequency response functions for a multi-degree-of-freedom
        !! system that has a general damping matrix, and is not necessarily 
        !! symmetric.  The problem is treated as the solution to the linear
        !! system \( \left( K - \omega^{2} M + j \omega C \right) H(\omega) = 
        !! F(\omega) \).
        real(real64), intent(in), dimension(:,:) :: mass
            !! The N-by-N mass matrix.
        real(real64), intent(in), dimension(:,:) :: damp
            !! The N-by-N damping matrix.
        real(real64), intent(in), dimension(:,:) :: stiff
            !! The N-by-N stiffness matrix.
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
        integer(int32), intent(out), optional, dimension(:) :: ranks
            !! Provides information on the rank of the dynamic stiffness matrix
            !! for each frequency.  If provided, this array must be the same
            !! length as freq.
        class(*), intent(inout), optional :: args
            !! An optional argument that can be used to communicate with
            !! the outside world.
        type(frf) :: rst
            !! The resulting frequency responses.

        ! Local Variables
        integer(int32) :: i, flag
        real(real64) :: df
        real(real64), allocatable, dimension(:) :: freq

        ! Input Checking
        if (abs(freq1 - freq2) < sqrt(epsilon(freq1))) error stop DYN_INVALID_INPUT_ERROR
        if (nfreq < 2) error stop DYN_INVALID_INPUT_ERROR

        ! Process
        df = (freq2 - freq1) / (nfreq - 1.0d0)
        allocate(freq(nfreq))
        freq = (/ (df * i + freq1, i = 0, nfreq - 1) /)
        rst = frequency_response(mass, damp, stiff, freq, frc, ranks = ranks, &
            args = args)
    end function

! ------------------------------------------------------------------------------
end module
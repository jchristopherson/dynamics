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
module dynamics_modal_analysis
    use iso_fortran_env
    use dynamics_error_handling
    use dynamics_helper
    use linalg, only : eigen, sort, csr_matrix, eigen, size
    implicit none
    private
    public :: compute_modal_damping
    public :: modal_response
    public :: normalize_mode_shapes

    interface modal_response
        module procedure :: modal_response_dense
        module procedure :: modal_response_sparse
    end interface

contains

! ------------------------------------------------------------------------------
    pure elemental function compute_modal_damping(lambda, alpha, beta) &
        result(rst)
        !! Computes the modal damping factors \( \zeta_i \) given the
        !! proportional damping terms \( \alpha \) and \( \beta \) where
        !! \( \alpha + \beta \omega_{i}^2 = 2 \zeta_{i} \omega_{i} \),
        !! \( \lambda_{i} = \omega_{i}^2 \), and \( \lambda_i \) is the
        !! \( i^{th} \) eigenvalue of the system.
        !! Equivalently, \(\zeta_i=(\alpha+\beta\omega_i^2)/(2\omega_i)\).
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
    pure subroutine modal_response_dense(mass, stiff, freqs, modeshapes)
        !! Computes the modal frequencies and modes shapes for 
        !! multi-degree-of-freedom system.
        !! The generalized eigenproblem is
        !! $$ K\boldsymbol{\phi}_i=\lambda_iM\boldsymbol{\phi}_i,
        !! \qquad \omega_i=\sqrt{\lambda_i}. $$
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

        ! Local Variables
        integer(int32) :: n
        complex(real64), allocatable, dimension(:) :: vals
        complex(real64), allocatable, dimension(:,:) :: vecs
        
        ! Initialization
        n = size(mass, 1)

        ! Input Checking
        if (n < 1) error stop DYN_INVALID_INPUT_ERROR
        if (size(mass, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR
        if (size(stiff, 1) /= size(stiff, 2)) error stop DYN_MATRIX_SIZE_ERROR
        if (size(stiff, 1) /= n .or. size(stiff, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR
        if (.not.is_symmetric(mass) .or. .not.is_symmetric(stiff)) &
            error stop DYN_INVALID_INPUT_ERROR

        ! Memory allocations
        allocate(vals(n))
        if (present(modeshapes)) allocate(vecs(n, n))

        ! Solve the eigen problem
        if (present(modeshapes)) then
            call eigen(stiff, mass, vals, rvecs = vecs)
            call sort(vals, vecs)
            allocate(modeshapes(n, n), source = real(vecs))
        else
            call eigen(stiff, mass, vals)
            call sort(vals)
        end if

        ! Convert the eigenvalues to frequency values
        if (any(real(vals) <= 0.0d0)) error stop DYN_INVALID_INPUT_ERROR
        allocate(freqs(n), source = sqrt(real(vals)))
    end subroutine

! ------------------------------------------------------------------------------
    subroutine modal_response_sparse(mass, stiff, nmodes, freqs, modeshapes)
        !! Computes selected modal frequencies and mode shapes for a
        !! multi-degree-of-freedom system using CSR sparse matrices.
        !! The generalized eigenproblem is
        !! $$ K\boldsymbol{\phi}_i=\lambda_iM\boldsymbol{\phi}_i,
        !! \qquad \omega_i=\sqrt{\lambda_i}. $$
        type(csr_matrix), intent(in) :: mass
            !! The N-by-N symmetric positive-definite mass matrix.
        type(csr_matrix), intent(in) :: stiff
            !! The N-by-N symmetric stiffness matrix.
        integer(int32), intent(in) :: nmodes
            !! The number of lowest-frequency modes to compute.  This value
            !! must be greater than zero and less than N.
        real(real64), intent(out), allocatable, dimension(:) :: freqs
            !! An allocatable NMODES-element array containing the modal
            !! frequencies in ascending order with units of rad/s.
        real(real64), intent(out), allocatable, optional, dimension(:,:) :: &
            modeshapes
            !! An optional, allocatable N-by-NMODES matrix containing one mode
            !! shape per column.

        integer(int32) :: i, j, loc, n
        real(real64) :: mass_scale, stiff_scale, temp, tol
        real(real64), allocatable, dimension(:) :: temp_vec, vals
        real(real64), allocatable, dimension(:,:) :: vecs

        n = size(mass, 1)

        if (n < 1) error stop DYN_INVALID_INPUT_ERROR
        if (size(mass, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR
        if (size(stiff, 1) /= n .or. size(stiff, 2) /= n) &
            error stop DYN_MATRIX_SIZE_ERROR
        if (nmodes < 1 .or. nmodes >= n) error stop DYN_INVALID_INPUT_ERROR

        tol = 10.0d0 * epsilon(0.0d0)
        mass_scale = 1.0d0
        stiff_scale = 1.0d0
        if (size(mass%values) > 0) &
            mass_scale = max(mass_scale, maxval(abs(mass%values)))
        if (size(stiff%values) > 0) &
            stiff_scale = max(stiff_scale, maxval(abs(stiff%values)))
        do j = 1, n
            do i = mass%row_indices(j), mass%row_indices(j + 1) - 1
                if (abs(mass%values(i) - &
                    mass%get(mass%column_indices(i), j)) > &
                    tol * mass_scale) error stop DYN_INVALID_INPUT_ERROR
            end do
            do i = stiff%row_indices(j), stiff%row_indices(j + 1) - 1
                if (abs(stiff%values(i) - &
                    stiff%get(stiff%column_indices(i), j)) > &
                    tol * stiff_scale) error stop DYN_INVALID_INPUT_ERROR
            end do
        end do

        if (present(modeshapes)) then
            call eigen(stiff, mass, nmodes, vals, vecs, sigma = 0.0d0)
            allocate(temp_vec(n))
            do i = 1, size(vals) - 1
                loc = i - 1 + minloc(vals(i:), 1)
                if (loc /= i) then
                    temp = vals(i)
                    vals(i) = vals(loc)
                    vals(loc) = temp
                    temp_vec = vecs(:,i)
                    vecs(:,i) = vecs(:,loc)
                    vecs(:,loc) = temp_vec
                end if
            end do
            allocate(modeshapes(n, size(vals)), source = vecs)
        else
            call eigen(stiff, mass, nmodes, vals, sigma = 0.0d0)
            call sort(vals)
        end if

        if (size(vals) /= nmodes .or. any(vals <= 0.0d0)) &
            error stop DYN_INVALID_INPUT_ERROR
        allocate(freqs(nmodes), source = sqrt(vals))
    end subroutine

! ------------------------------------------------------------------------------
    pure subroutine normalize_mode_shapes(x)
        !! Normalizes mode shape vectors such that the largest magnitude
        !! value in the vector is one.
        !! For each column \(\boldsymbol{\phi}_i\), the operation is
        !! $$ \boldsymbol{\phi}_i\leftarrow
        !! \frac{\boldsymbol{\phi}_i}{\phi_{i,k}},\qquad
        !! k=\arg\max_j|\phi_{i,j}|. $$
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

! ------------------------------------------------------------------------------
end module
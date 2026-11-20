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
module dynamics_flexible_body_solvers
    !! Nonlinear static equilibrium for flexible finite-element meshes.
    !!
    !! The solver uses the mesh's internal force and consistent tangent, so
    !! geometric and material nonlinearity are those implemented by its
    !! elements and materials. External forces are fixed in the global frame.
    !! This module does not integrate in time, follow equilibrium branches,
    !! or manage history-dependent material state.
    use iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use dynamics_error_handling, only : DYN_ARRAY_SIZE_ERROR, DYN_INVALID_INPUT_ERROR, &
        DYN_CONVERGENCE_ERROR, DYN_NONMONOTONIC_ARRAY_ERROR
    use dynamics_flexible_assembly, only : finite_element_mesh
    use dynamics_boundary_conditions, only : nodal_condition, resolve_boundary_conditions
    use dynamics_structural, only : apply_boundary_conditions
    use linalg, only : csr_matrix, msr_matrix, lu_factor, pgmres_solve, matmul
    use lapack, only : DGETRF, DGETRS
    implicit none
    private
    public :: flexible_static_settings
    public :: flexible_static_info
    public :: solve_flexible_static

    type :: flexible_static_settings
        !! Controls the residual test and Newton iteration. Default values are
        !! used when [[solve_flexible_static]] is called without settings.
        logical :: sparse = .false.
            !! If true, assemble a CSR tangent and solve with ILUT-preconditioned
            !! GMRES. Otherwise use dense LU.
        integer :: max_iterations = 30
            !! Maximum number of Newton corrections; must be positive.
        integer :: max_backtracks = 20
            !! Maximum halvings of the correction after the full step is
            !! tried; zero permits only a full Newton step.
        integer :: nthreads = 1
            !! Positive number of threads requested for force and tangent
            !! assembly. Parallel evaluation requires thread-safe elements
            !! and materials, as documented by [[finite_element_mesh]].
        real(real64) :: absolute_tolerance = 1.0d-10
            !! Nonnegative absolute free-residual norm tolerance. This and
            !! `relative_tolerance` cannot both be zero.
        real(real64) :: relative_tolerance = 1.0d-8
            !! Nonnegative tolerance relative to the larger of the initial
            !! free-residual norm and the applied free-force norm.
    end type

    type :: flexible_static_info
        !! Diagnostics for the initial or last accepted iterate, even when
        !! the solve fails and the output displacement is left unchanged.
        logical :: converged = .false.
            !! True only when the residual test passes or all DOFs are fixed.
        logical :: linear_solve_failed = .false.
            !! True if the tangent, factorization, or correction solve failed
            !! its validity or true linear-residual check. False does not
            !! imply convergence: line search or the iteration limit can fail.
        integer :: iterations = 0
            !! Number of accepted Newton corrections, not force evaluations
            !! or backtracking trials.
        real(real64) :: residual_norm = 0.0d0
            !! Euclidean norm of `force - internal_force` on free DOFs at the
            !! last accepted iterate. It can be nonfinite if evaluation fails;
            !! it is zero when all DOFs are prescribed.
    end type

contains
    subroutine solve_flexible_static(mesh, force, conditions, displacement, settings, info)
        !! Solve static equilibrium on the free equations of a flexible mesh:
        !!
        !! $$r_f(q) = f_f - f_{\mathrm{int},f}(q) = 0.$$
        !!
        !! Prescribed displacement values are imposed before evaluating the
        !! residual. At each iteration the assembled free tangent solves
        !! \(K_{ff}\Delta q_f = r_f\). A residual-norm backtracking search
        !! tries the full correction followed by successive halvings. A trial
        !! is accepted when it satisfies the convergence tolerance or reduces
        !! the residual norm by at least `1.0d-4 * alpha` of its previous value.
        !!
        !! The free residual norm must be at most
        !! `absolute_tolerance + relative_tolerance * max(initial_residual_norm, norm2(force_free))`.
        !! The norm is unscaled across DOFs; force and moment components may
        !! need a nondimensionalized model. If all DOFs are prescribed, the
        !! solve succeeds without evaluating forces or tangents.
        !!
        !! Element internal forces must use the positive internal-force
        !! convention and their tangents must be consistent derivatives.
        !! Evaluations must be safe to repeat at trial displacements; this
        !! routine does not commit or roll back history-dependent materials.
        !! Force is independent of displacement: follower loads and their
        !! tangents are not supported. For prescribed-load increments, call
        !! the routine repeatedly with the previous solution as the next guess.
        !!
        !! With `info`, nonlinear or linear-solver failure is reported there
        !! and `displacement` remains unchanged. Without `info`, such failure
        !! causes `error stop DYN_CONVERGENCE_ERROR`. Invalid inputs and errors
        !! raised by the mesh or elements cause an error stop regardless of
        !! whether `info` is supplied.
        type(finite_element_mesh), intent(inout) :: mesh
            !! Mesh to assemble. Its element workspaces may be modified; do
            !! not assemble or solve the same mesh concurrently from callers.
        real(real64), intent(in), dimension(:) :: force
            !! Finite, displacement-independent external force in full mesh
            !! equation order. Loads on prescribed DOFs are ignored by the
            !! equilibrium test but can be used when computing reactions.
        type(nodal_condition), intent(in), dimension(:) :: conditions
            !! Prescribed displacements addressed by node-table position and
            !! physical DOF label. Values must be finite and DOFs distinct.
            !! An empty array imposes no displacement constraints.
        real(real64), intent(inout), dimension(:) :: displacement
            !! Finite full displacement vector in mesh equation order.
            !! Its free entries are the initial guess; prescribed entries are
            !! replaced by condition values on success. On recoverable failure
            !! the entire input vector remains unchanged.
        type(flexible_static_settings), intent(in), optional :: settings
            !! Optional iteration and linear-solver settings. Omit to use
            !! [[flexible_static_settings]] defaults.
        type(flexible_static_info), intent(out), optional :: info
            !! Optional convergence diagnostics. Supplying this argument makes
            !! nonlinear and linear-solver failures recoverable.

        type(flexible_static_settings) :: opts
        type(flexible_static_info) :: status
        type(csr_matrix) :: sparse_k
        integer, allocatable, dimension(:) :: constrained
        integer, allocatable, dimension(:) :: sorted
        integer, allocatable, dimension(:) :: free
        integer, allocatable, dimension(:) :: pivots
        logical, allocatable, dimension(:) :: is_free
        real(real64), allocatable, dimension(:) :: q
        real(real64), allocatable, dimension(:) :: trial
        real(real64), allocatable, dimension(:) :: internal
        real(real64), allocatable, dimension(:) :: residual
        real(real64), allocatable, dimension(:) :: delta
        real(real64), allocatable, dimension(:,:) :: dense_k
        real(real64), allocatable, dimension(:,:) :: original_k
        real(real64) :: initial_norm
        real(real64) :: threshold
        real(real64) :: trial_norm
        real(real64) :: alpha
        integer :: ndof
        integer :: nfree
        integer :: i
        integer :: iteration
        integer :: backtrack
        integer :: lapack_info
        logical :: linear_ok
        logical :: accepted

        if (present(settings)) opts = settings
        if (opts%max_iterations < 1 .or. opts%max_backtracks < 0 .or. opts%nthreads < 1 .or. &
            .not.ieee_is_finite(opts%absolute_tolerance) .or. &
            .not.ieee_is_finite(opts%relative_tolerance)) error stop DYN_INVALID_INPUT_ERROR
        if (opts%absolute_tolerance < 0.0d0 .or. opts%relative_tolerance < 0.0d0 .or. &
            (opts%absolute_tolerance == 0.0d0 .and. opts%relative_tolerance == 0.0d0)) &
            error stop DYN_INVALID_INPUT_ERROR
        call resolve_boundary_conditions(mesh, conditions, constrained, ndof)
        if (size(force) /= ndof .or. size(displacement) /= ndof) error stop DYN_ARRAY_SIZE_ERROR
        if (ndof < 1 .or. any(.not.ieee_is_finite(force)) .or. &
            any(.not.ieee_is_finite(displacement)) .or. &
            any(.not.ieee_is_finite(conditions%value))) error stop DYN_INVALID_INPUT_ERROR

        ! Keep prescribed values separate from the free Newton corrections.
        allocate(is_free(ndof), source=.true.)
        do i = 1, size(constrained)
            if (.not.is_free(constrained(i))) error stop DYN_NONMONOTONIC_ARRAY_ERROR
            is_free(constrained(i)) = .false.
        end do
        free = pack([(i, i=1,ndof)], is_free)
        nfree = size(free)

        ! Work on a private candidate so failure never publishes a partial step.
        q = displacement
        q(constrained) = conditions%value
        if (nfree == 0) then
            status%converged = .true.
            displacement = q
            if (present(info)) info = status
            return
        end if
        call mesh%assemble_internal_force(q, internal, nthreads=opts%nthreads)
        residual = force(free) - internal(free)
        status%residual_norm = norm2(residual)
        if (.not.all(ieee_is_finite(residual))) then
            call failure()
            return
        end if
        initial_norm = status%residual_norm

        ! Use a fixed reference throughout this solve, including displacement-driven cases.
        threshold = opts%absolute_tolerance + opts%relative_tolerance * max(initial_norm, norm2(force(free)))
        if (.not.ieee_is_finite(threshold) .or. .not.ieee_is_finite(initial_norm)) then
            call failure()
            return
        end if

        do iteration = 1, opts%max_iterations
            if (status%residual_norm <= threshold) exit
            if (opts%sparse) then
                call mesh%assemble_stiffness(q, sparse_k, nthreads=opts%nthreads)
                if (size(constrained) > 0) then
                    ! Prescribed values already enter q; their Newton corrections are zero.
                    sorted = constrained
                    sparse_k = apply_boundary_conditions(sorted, sparse_k)
                end if
                call sparse_correction(sparse_k, residual, delta, linear_ok)
            else
                call mesh%assemble_stiffness(q, dense_k, nthreads=opts%nthreads)
                dense_k = dense_k(free,free)
                linear_ok = all(ieee_is_finite(dense_k))
                if (linear_ok) then
                    ! LU overwrites the tangent; retain it to check the true residual.
                    original_k = dense_k
                    allocate(pivots(nfree))
                    call DGETRF(nfree, nfree, dense_k, nfree, pivots, lapack_info)
                    linear_ok = lapack_info == 0
                    if (linear_ok) then
                        delta = residual
                        call DGETRS('N', nfree, 1, dense_k, nfree, pivots, delta, nfree, lapack_info)
                        linear_ok = lapack_info == 0 .and. all(ieee_is_finite(delta))
                        if (linear_ok) linear_ok = norm2(residual - matmul(original_k, delta)) <= &
                            1.0d-10 * norm2(residual)
                    end if
                    deallocate(pivots)
                end if
            end if
            if (.not.linear_ok) then
                status%linear_solve_failed = .true.
                call failure()
                return
            end if
            accepted = .false.
            alpha = 1.0d0

            ! Test equilibrium at the full trial before halving the correction.
            do backtrack = 0, opts%max_backtracks
                trial = q
                trial(free) = q(free) + alpha * delta
                if (all(ieee_is_finite(trial))) then
                    call mesh%assemble_internal_force(trial, internal, nthreads=opts%nthreads)
                    if (all(ieee_is_finite(internal(free)))) then
                        trial_norm = norm2(force(free) - internal(free))
                        if (ieee_is_finite(trial_norm)) then
                            if (trial_norm <= threshold .or. &
                                trial_norm <= (1.0d0 - 1.0d-4 * alpha) * status%residual_norm) then
                                accepted = .true.
                                exit
                            end if
                        end if
                    end if
                end if
                alpha = 0.5d0 * alpha
            end do
            if (.not.accepted) then
                call failure()
                return
            end if
            q = trial
            residual = force(free) - internal(free)
            status%residual_norm = trial_norm
            status%iterations = iteration
        end do
        if (status%residual_norm > threshold) then
            call failure()
            return
        end if
        status%converged = .true.
        displacement = q
        if (present(info)) info = status
    contains
        subroutine failure()
            !! Report a recoverable failure or terminate when no status was requested.
            if (present(info)) then
                info = status
            else
                error stop DYN_CONVERGENCE_ERROR
            end if
        end subroutine
    end subroutine

    subroutine sparse_correction(k, rhs, delta, success)
        !! Solve a reduced CSR Newton system using ILUT-preconditioned
        !! restarted GMRES and verify the true residual to relative tolerance
        !! 1e-10. ILUT uses a drop tolerance of sqrt(machine epsilon).
        type(csr_matrix), intent(in) :: k
            !! Nonempty reduced tangent matrix.
        real(real64), intent(in), dimension(:) :: rhs
            !! Free residual vector defining the correction right-hand side.
        real(real64), allocatable, intent(out), dimension(:) :: delta
            !! Computed correction on success; not usable when `success` is
            !! false, even if allocated.
        logical, intent(out) :: success
            !! True only if GMRES and the true-residual check both succeed.

        type(msr_matrix) :: preconditioner
        integer, allocatable, dimension(:) :: upper_start
        integer :: n
        integer :: flag

        success = .false.
        n = size(rhs)
        if (any(.not.ieee_is_finite(k%values))) return
        allocate(upper_start(n))
        
        ! The factors precondition this tangent; GMRES still solves the full system.
        call lu_factor(k, preconditioner, upper_start, droptol=sqrt(epsilon(1.0d0)), flag=flag)
        if (flag /= 0) return
        if (any(.not.ieee_is_finite(preconditioner%values))) return
        call pgmres_solve(k, preconditioner, upper_start, rhs, delta, &
            im=min(n, 50), tol=1.0d-10, maxits=500, flag=flag)
        if (flag /= 0 .or. .not.allocated(delta)) return
        if (any(.not.ieee_is_finite(delta))) return
        success = norm2(rhs - matmul(k, delta)) <= 1.0d-10 * norm2(rhs)
    end subroutine
end module

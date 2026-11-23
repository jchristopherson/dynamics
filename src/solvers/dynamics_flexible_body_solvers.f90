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
    !! Static equilibrium and transient dynamics for flexible finite-element meshes.
    !!
    !! Both solvers use the mesh's internal force. Static equilibrium also
    !! uses the assembled tangent; transient integration uses a constant
    !! reference mass and lets the integrator approximate its own Jacobian.
    !! Neither solver manages history-dependent material state.
    use iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use dynamics_error_handling, only : DYN_ARRAY_SIZE_ERROR, DYN_INVALID_INPUT_ERROR, &
        DYN_CONVERGENCE_ERROR, DYN_NONMONOTONIC_ARRAY_ERROR, DYN_CONSTRAINT_ERROR
    use dynamics_flexible_assembly, only : finite_element_mesh
    use dynamics_boundary_conditions, only : nodal_condition, resolve_boundary_conditions
    use dynamics_structural, only : apply_boundary_conditions
    use dynamics_solver_reporting, only : static_solver_progress, static_solver_callback
    use linalg, only : csr_matrix, msr_matrix, lu_factor, pgmres_solve, matmul
    use lapack, only : DGETRF, DGETRS, DPOTRF, DGESVD
    use diffeq, only : ode_container, ggl_container, ode_integrator, bdf, rosenbrock, &
        kennedy_carpenter_4, kennedy_carpenter_5
    implicit none
    private
    public :: flexible_static_settings
    public :: flexible_static_info
    public :: solve_flexible_static
    public :: solve_flexible_transient, flexible_transient_settings, flexible_load_callback
    public :: FLEXIBLE_BDF, FLEXIBLE_ROSENBROCK, FLEXIBLE_KENNEDY_CARPENTER_4, FLEXIBLE_KENNEDY_CARPENTER_5

    integer, parameter :: FLEXIBLE_BDF = 1
        !! Backward differentiation formula (the transient default).
    integer, parameter :: FLEXIBLE_ROSENBROCK = 2
        !! Linearly implicit Rosenbrock method.
    integer, parameter :: FLEXIBLE_KENNEDY_CARPENTER_4 = 3
        !! Fourth-order Kennedy-Carpenter implicit Runge-Kutta method.
    integer, parameter :: FLEXIBLE_KENNEDY_CARPENTER_5 = 4
        !! Fifth-order Kennedy-Carpenter implicit Runge-Kutta method.

    type :: flexible_transient_settings
        !! Controls for [[solve_flexible_transient]]. All real-valued
        !! settings must be finite.
        integer :: solver = FLEXIBLE_BDF
            !! One of `FLEXIBLE_BDF`, `FLEXIBLE_ROSENBROCK`,
            !! `FLEXIBLE_KENNEDY_CARPENTER_4`, or `FLEXIBLE_KENNEDY_CARPENTER_5`.
        real(real64) :: absolute_tolerance = 1.0d-8
            !! Positive absolute error tolerance applied to every integrated
            !! state component, including algebraic variables when present.
        real(real64) :: relative_tolerance = 1.0d-6
            !! Nonnegative relative error tolerance applied to every state
            !! component. Models with disparate scales may need rescaling.
        real(real64) :: minimum_step_size = 10.0d0 * epsilon(1.0d0)
            !! Positive lower step-size bound in the units of `times`.
        real(real64) :: maximum_step_size = huge(1.0d0)
            !! Upper step-size bound, at least `minimum_step_size`. It does
            !! not set the spacing of the requested output samples.
        integer :: step_limit = 1000000
            !! Positive integration step limit passed to the selected solver.
        integer :: nthreads = 1
            !! Positive thread count requested for internal-force assembly.
            !! Parallel assembly requires thread-safe elements and materials.
        real(real64) :: constraint_tolerance = 1.0d-8
            !! Nonnegative absolute bound on each initial entry of `g` and
            !! `G*v`. This is not an integration or drift tolerance.
    end type

    abstract interface
        subroutine flexible_load_callback(time, force, args)
            !! Supply an applied load that depends on time, not the state.
            !! Evaluations may repeat or be made out of chronological order.
            !! Smooth loads are required; do not advance external state here.
            import :: real64
            real(real64), intent(in) :: time
                !! Internal evaluation time, not necessarily an output time.
            real(real64), intent(out), dimension(:) :: force
                !! Fill every entry with a finite value in full mesh equation
                !! order, including zeros for unloaded DOFs. Prescribed DOF
                !! entries do not enter the free equations.
            class(*), intent(in), optional :: args
                !! User data forwarded only if supplied to the solve.
        end subroutine
    end interface

    type :: flexible_transient_context
        !! Per-solve data shared by the private integration callbacks.
        type(finite_element_mesh), pointer :: mesh => null()
        procedure(flexible_load_callback), pointer, nopass :: load => null()
        class(*), pointer :: load_args => null()
        integer, allocatable, dimension(:) :: free
            !! Maps reduced coordinates back to full mesh equations.
        real(real64), allocatable, dimension(:) :: support
            !! Full displacement template with zero entries at free DOFs.
        real(real64), allocatable, dimension(:,:) :: mass
            !! Constant free-free block of the assembled reference mass.
        integer :: nthreads
    end type

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
    subroutine solve_flexible_transient(mesh, times, conditions, initial_displacement, initial_velocity, &
        displacement, velocity, settings, load, args, multipliers)
        !! Integrate a mesh with a dense, constant reference mass. Prescribed
        !! displacements are constant and eliminated from the integration
        !! equations. The reduced mass must be positive definite. If all
        !! DOFs are prescribed and there are no joints, return the prescribed
        !! history without assembling mass or evaluating loads.
        !!
        !! For meshes without joints the free equations are `qdot = v` and
        !! `M*vdot = f_ext - f_int - f_damp`. For meshes with joints a private stabilized
        !! GGL adapter adds the joint force `-transpose(G)*lambda`, where `G`
        !! is the constraint Jacobian restricted to free columns. It also
        !! adds a kinematic correction to control drift.
        !!
        !! With joints, the initial state after imposing supports must satisfy
        !! both `g(q) = 0` and `G(q)*v = 0` within `constraint_tolerance`.
        !! The initial reduced Jacobian must have independent rows. The
        !! solver computes consistent multipliers but does not project the
        !! input displacement or velocity. Avoid configurations where the
        !! joint Jacobian loses rank during the motion.
        !!
        !! Invalid inputs and integration failures cause `error stop`.
        !! Moving supports, state-dependent applied loads and inertia, and
        !! history-dependent materials are not supported.
        type(finite_element_mesh), intent(inout), target :: mesh
            !! Mesh to assemble. Element workspaces may change; do not
            !! assemble or solve this mesh concurrently from callers.
        real(real64), intent(in), dimension(:) :: times
            !! Finite, strictly increasing output times, with at least two
            !! entries. Integration starts at `times(1)`.
        type(nodal_condition), intent(in), dimension(:) :: conditions
            !! Constant prescribed displacements. Values must be finite and
            !! DOFs distinct; an empty array imposes no supports.
        real(real64), intent(in), dimension(:) :: initial_displacement
            !! Finite initial displacement in full mesh equation order.
            !! Prescribed entries are ignored in favor of `conditions`.
        real(real64), intent(in), dimension(:) :: initial_velocity
            !! Finite initial coordinate rates in full mesh equation order.
            !! Prescribed entries are ignored; those rates are set to zero.
        real(real64), allocatable, intent(out), dimension(:,:) :: displacement
            !! Full displacement history, indexed by (mesh DOF, time sample).
        real(real64), allocatable, intent(out), dimension(:,:) :: velocity
            !! Full coordinate-rate history, indexed by (mesh DOF, time sample).
            !! Rotational entries are rotation-vector rates, not angular velocity.
        type(flexible_transient_settings), intent(in), optional :: settings
            !! Omit to use [[flexible_transient_settings]] defaults.
        procedure(flexible_load_callback), optional :: load
            !! Smooth applied load. Omit to use zero applied force.
        class(*), intent(in), optional, target :: args
            !! Optional user data forwarded unchanged to `load`.
        real(real64), allocatable, intent(out), optional, dimension(:,:) :: multipliers
            !! Physical joint force multipliers `lambda` in constraint assembly
            !! row order, with one column per time sample. Allocated with zero
            !! rows if the mesh has no joints. The kinematic correction
            !! multiplier is not returned.

        type(flexible_transient_settings) :: opts
        type(flexible_transient_context) :: ctx
        class(ode_integrator), allocatable :: solver
        type(ode_container) :: ode
        type(ggl_container) :: dae
        integer, allocatable, dimension(:) :: constrained
        logical, allocatable, dimension(:) :: is_free
        real(real64), allocatable, dimension(:,:) :: full_mass
        real(real64), allocatable, dimension(:,:) :: factor
        real(real64), allocatable, dimension(:,:) :: sol
        real(real64), allocatable, dimension(:,:) :: gradient
        real(real64), allocatable, dimension(:) :: state
        integer :: ndof
        integer :: n
        integer :: nc
        integer :: nt
        integer :: i
        integer :: row
        integer :: flag

        if (present(settings)) opts = settings
        if (.not.all(ieee_is_finite([opts%absolute_tolerance, opts%relative_tolerance, opts%minimum_step_size, &
            opts%maximum_step_size, opts%constraint_tolerance]))) error stop DYN_INVALID_INPUT_ERROR
        if (opts%absolute_tolerance <= 0.0d0 .or. opts%relative_tolerance < 0.0d0 .or. &
            opts%minimum_step_size <= 0.0d0 .or. opts%maximum_step_size < opts%minimum_step_size .or. &
            opts%constraint_tolerance < 0.0d0 .or. opts%step_limit < 1 .or. opts%nthreads < 1) &
            error stop DYN_INVALID_INPUT_ERROR
        select case (opts%solver)
        case (FLEXIBLE_BDF)
            allocate(bdf :: solver)
        case (FLEXIBLE_ROSENBROCK)
            allocate(rosenbrock :: solver)
        case (FLEXIBLE_KENNEDY_CARPENTER_4)
            allocate(kennedy_carpenter_4 :: solver)
        case (FLEXIBLE_KENNEDY_CARPENTER_5)
            allocate(kennedy_carpenter_5 :: solver)
        case default
            error stop DYN_INVALID_INPUT_ERROR
        end select

        nt = size(times)
        if (nt < 2) error stop DYN_ARRAY_SIZE_ERROR
        if (.not.all(ieee_is_finite(times))) error stop DYN_INVALID_INPUT_ERROR
        if (any(times(2:) <= times(:nt-1))) error stop DYN_NONMONOTONIC_ARRAY_ERROR
        call resolve_boundary_conditions(mesh, conditions, constrained, ndof)
        if (size(initial_displacement) /= ndof .or. size(initial_velocity) /= ndof) error stop DYN_ARRAY_SIZE_ERROR
        if (ndof < 1 .or. .not.all(ieee_is_finite(initial_displacement)) .or. &
            .not.all(ieee_is_finite(initial_velocity)) .or. .not.all(ieee_is_finite(conditions%value))) &
            error stop DYN_INVALID_INPUT_ERROR
        allocate(is_free(ndof), source=.true.)
        do i = 1, size(constrained)
            if (.not.is_free(constrained(i))) error stop DYN_NONMONOTONIC_ARRAY_ERROR
            is_free(constrained(i)) = .false.
        end do
        ctx%mesh => mesh
        if (present(load)) ctx%load => load
        if (present(args)) ctx%load_args => args
        ctx%nthreads = opts%nthreads
        ctx%free = pack([(i, i=1,ndof)], is_free)
        allocate(ctx%support(ndof), source=0.0d0)
        ctx%support(constrained) = conditions%value
        n = size(ctx%free)
        nc = mesh%number_constraints()
        if (nc > n) error stop DYN_CONSTRAINT_ERROR

        allocate(displacement(ndof,nt), velocity(ndof,nt), source=0.0d0)
        displacement = spread(ctx%support, 2, nt)
        if (present(multipliers)) allocate(multipliers(nc,nt), source=0.0d0)
        if (n == 0) return

        call mesh%assemble_mass(full_mass)
        ctx%mass = full_mass(ctx%free,ctx%free)
        deallocate(full_mass)
        if (.not.all(ieee_is_finite(ctx%mass))) error stop DYN_INVALID_INPUT_ERROR
        ! Validate definiteness without overwriting the mass used by callbacks.
        factor = ctx%mass
        call DPOTRF('L', n, factor, n, flag)
        if (flag /= 0) error stop DYN_INVALID_INPUT_ERROR
        deallocate(factor)

        call solver%set_absolute_tolerance(opts%absolute_tolerance)
        call solver%set_relative_tolerance(opts%relative_tolerance)
        call solver%set_minimum_step_size(opts%minimum_step_size)
        call solver%set_maximum_step_size(opts%maximum_step_size)
        call solver%set_step_limit(opts%step_limit)
        call solver%set_allow_overshoot(.false.)

        ! State order is [q, v] or [q, v, lambda, mu] with joints.
        allocate(state(2*n+2*nc), source=0.0d0)
        state(:n) = initial_displacement(ctx%free)
        state(n+1:2*n) = initial_velocity(ctx%free)
        if (nc == 0) then
            ode%fcn => transient_rhs
            ode%mass_matrix => transient_state_mass
            call ode%set_is_mass_matrix_dependent(.false.)
            call solver%solve(ode, times, state, args=ctx)
        else
            dae = ggl_container(n, nc)
            dae%fcn => transient_force
            dae%mass_matrix => transient_physical_mass
            dae%constraints => transient_constraints
            dae%constraint_jacobian => transient_gradient
            dae%constraint_curvature => transient_curvature
            call dae%set_is_mass_matrix_dependent(.false.)
            if (.not.dae%validate_initial_conditions(times(1), state, opts%constraint_tolerance, ctx)) &
                error stop DYN_CONSTRAINT_ERROR
            allocate(gradient(nc,n))
            call transient_gradient(times(1), state(:n), gradient, ctx)
            if (.not.independent_constraints(gradient)) error stop DYN_CONSTRAINT_ERROR
            call dae%initialize_multipliers(times(1), state, ctx)
            call solver%solve(dae, times, state, args=ctx)
        end if

        sol = solver%get_solution()
        if (size(sol,1) < nt .or. size(sol,2) /= size(state)+1) error stop DYN_CONVERGENCE_ERROR
        do i = 1, nt
            ! diffeq returns internal steps when given just the two endpoints.
            row = i
            if (nt == 2 .and. i == 2) row = size(sol,1)
            if (abs(sol(row,1) - times(i)) > spacing(max(abs(times(i)), abs(sol(row,1))))) &
                error stop DYN_CONVERGENCE_ERROR
            if (.not.all(ieee_is_finite(sol(row,:)))) error stop DYN_CONVERGENCE_ERROR
            displacement(ctx%free,i) = sol(row,2:n+1)
            velocity(ctx%free,i) = sol(row,n+2:2*n+1)
            if (present(multipliers)) multipliers(:,i) = sol(row,2*n+2:2*n+nc+1)
        end do
    end subroutine

    subroutine transient_force(t, y, f, args)
        !! Assemble the free net force from `[q, v]`. Joint
        !! forces are excluded here because the GGL adapter supplies them.
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: y
        real(real64), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        real(real64), allocatable, dimension(:) :: q
        real(real64), allocatable, dimension(:) :: v
        real(real64), allocatable, dimension(:) :: internal
        real(real64), allocatable, dimension(:) :: applied

        if (.not.present(args)) error stop DYN_INVALID_INPUT_ERROR
        select type (args)
        type is (flexible_transient_context)
            q = args%support
            q(args%free) = y(:size(f))
            allocate(v(size(q)), source=0.0d0)
            v(args%free) = y(size(f)+1:)
            call args%mesh%assemble_internal_force(q, internal, nthreads=args%nthreads, velocity=v)
            f = -internal(args%free)
            if (associated(args%load)) then
                allocate(applied(size(q)))
                if (associated(args%load_args)) then
                    call args%load(t, applied, args%load_args)
                else
                    call args%load(t, applied)
                end if
                if (.not.all(ieee_is_finite(applied))) error stop DYN_INVALID_INPUT_ERROR
                f = f + applied(args%free)
            end if
        class default
            error stop DYN_INVALID_INPUT_ERROR
        end select
    end subroutine

    subroutine transient_rhs(t, y, f, args)
        !! Unconstrained RHS for the mass-matrix system with state `[q, v]`.
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: y
        real(real64), intent(out), dimension(:) :: f
        class(*), intent(in), optional :: args
        integer :: n

        n = size(y) / 2
        f(:n) = y(n+1:)
        call transient_force(t, y, f(n+1:), args = args)
    end subroutine

    subroutine transient_physical_mass(t, y, mass, args)
        !! Return the cached reduced physical mass, independent of time and state.
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: y
        real(real64), intent(out), dimension(:,:) :: mass
        class(*), intent(in), optional :: args

        if (.not.present(args)) error stop DYN_INVALID_INPUT_ERROR
        select type (args)
        type is (flexible_transient_context)
            mass = args%mass
        class default
            error stop DYN_INVALID_INPUT_ERROR
        end select
    end subroutine

    subroutine transient_state_mass(t, y, mass, args)
        !! Unconstrained state mass `diag(I, M)`: unit kinematic block and
        !! the reduced physical mass in the momentum block.
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: y
        real(real64), intent(out), dimension(:,:) :: mass
        class(*), intent(in), optional :: args
        integer :: n
        integer :: i

        n = size(y) / 2
        mass = 0.0d0
        do i = 1, n
            mass(i,i) = 1.0d0
        end do
        call transient_physical_mass(t, y(:n), mass(n+1:,n+1:), args = args)
    end subroutine

    subroutine transient_joint_data(q, g, jac, args)
        !! Evaluate the joints in full coordinates, including the prescribed
        !! displacements, then retain only free columns in the Jacobian.
        real(real64), intent(in), dimension(:) :: q
        real(real64), allocatable, intent(out), dimension(:) :: g
        real(real64), allocatable, intent(out), dimension(:,:) :: jac
        class(*), intent(in), optional :: args
        real(real64), allocatable, dimension(:) :: full_q
        real(real64), allocatable, dimension(:,:) :: full_jac

        if (.not.present(args)) error stop DYN_INVALID_INPUT_ERROR
        select type (args)
        type is (flexible_transient_context)
            full_q = args%support
            full_q(args%free) = q
            call args%mesh%assemble_constraints(full_q, g, full_jac)
            jac = full_jac(:,args%free)
        class default
            error stop DYN_INVALID_INPUT_ERROR
        end select
    end subroutine

    subroutine transient_constraints(t, q, g, args)
        !! GGL callback for the autonomous joint position residual `g(q)`.
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:) :: g
        class(*), intent(in), optional :: args
        real(real64), allocatable, dimension(:) :: values
        real(real64), allocatable, dimension(:,:) :: jac

        call transient_joint_data(q, values, jac, args = args)
        g = values
    end subroutine

    subroutine transient_gradient(t, q, jac, args)
        !! GGL callback for `G = dg/dq` with respect to free coordinates.
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(out), dimension(:,:) :: jac
        class(*), intent(in), optional :: args
        real(real64), allocatable, dimension(:) :: values
        real(real64), allocatable, dimension(:,:) :: gradient

        call transient_joint_data(q, values, gradient, args = args)
        jac = gradient
    end subroutine

    subroutine transient_curvature(t, q, direction, v, curvature, args)
        !! Approximate `DG(q)[direction]*v` by centered differences of the
        !! reduced Jacobian. Normalize the direction before perturbing `q`
        !! so the displacement increment is set by the coordinate scale.
        real(real64), intent(in) :: t
        real(real64), intent(in), dimension(:) :: q
        real(real64), intent(in), dimension(:) :: direction
        real(real64), intent(in), dimension(:) :: v
        real(real64), intent(out), dimension(:) :: curvature
        class(*), intent(in), optional :: args
        real(real64), dimension(size(curvature),size(q)) :: plus
        real(real64), dimension(size(curvature),size(q)) :: minus
        real(real64) :: scale
        real(real64) :: h

        scale = maxval(abs(direction))
        curvature = 0.0d0
        if (scale == 0.0d0 .or. all(v == 0.0d0)) return
        ! Off the manifold, direction is qdot, not necessarily v.
        h = epsilon(1.0d0)**(1.0d0/3.0d0) * max(1.0d0, maxval(abs(q)))
        call transient_gradient(t, q + h * (direction / scale), plus, args = args)
        call transient_gradient(t, q - h * (direction / scale), minus, args = args)
        curvature = (scale / (2.0d0*h)) * matmul(plus - minus, v)
    end subroutine

    logical function independent_constraints(gradient) result(ok)
        !! Test initial reduced row rank by SVD. With no more rows than
        !! columns, the cutoff is `epsilon * size(gradient,2) * sigma_max`.
        real(real64), intent(in), dimension(:,:) :: gradient
        real(real64), dimension(size(gradient,1),size(gradient,2)) :: copy
        real(real64), dimension(size(gradient,1)) :: sigma
        real(real64), dimension(1,1) :: unused_u
        real(real64), dimension(1,1) :: unused_v
        real(real64), allocatable, dimension(:) :: work
        integer :: m
        integer :: n
        integer :: flag

        m = size(gradient,1)
        n = size(gradient,2)
        copy = gradient
        allocate(work(max(3*m+n,5*m)))
        call DGESVD('N', 'N', m, n, copy, m, sigma, unused_u, 1, unused_v, 1, work, size(work), flag)
        ok = .false.
        if (flag == 0) ok = sigma(m) > epsilon(1.0d0) * n * sigma(1)
    end function

    subroutine solve_flexible_static(mesh, force, conditions, displacement, settings, info, &
        progress, args, multipliers)
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
        !! Joints in `mesh%joints` are enforced with Lagrange multipliers
        !! \(\lambda\) in the augmented Newton system
        !! \([K_{ff}\ G_f^T; G_f\ 0]\). The tangent omits the multiplier-weighted
        !! constraint Hessian, so the iteration is quasi-Newton. Joints require
        !! the dense path and at least one free DOF.
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
        procedure(static_solver_callback), optional :: progress
            !! Optional callback invoked after each accepted Newton iterate.
        class(*), intent(in), optional :: args
            !! Optional user-supplied data forwarded to the progress callback.
        real(real64), allocatable, intent(out), optional, dimension(:) :: multipliers
            !! Optional converged Lagrange multipliers in joint-row order.
            !! Allocated only on success; empty when there are no joints.

        type(flexible_static_settings) :: opts
        type(flexible_static_info) :: status
        type(static_solver_progress) :: report
        type(csr_matrix) :: sparse_k
        integer, allocatable, dimension(:) :: constrained
        integer, allocatable, dimension(:) :: sorted
        integer, allocatable, dimension(:) :: free
        integer, allocatable, dimension(:) :: pivots
        logical, allocatable, dimension(:) :: is_free
        real(real64), allocatable, dimension(:) :: q
        real(real64), allocatable, dimension(:) :: trial
        real(real64), allocatable, dimension(:) :: residual
        real(real64), allocatable, dimension(:) :: trial_residual
        real(real64), allocatable, dimension(:) :: delta
        real(real64), allocatable, dimension(:) :: rhs
        real(real64), allocatable, dimension(:) :: lam
        real(real64), allocatable, dimension(:) :: lam_trial
        real(real64), allocatable, dimension(:) :: gcon
        real(real64), allocatable, dimension(:,:) :: dense_k
        real(real64), allocatable, dimension(:,:) :: jcon
        real(real64), allocatable, dimension(:,:) :: augmented
        real(real64), allocatable, dimension(:,:) :: original_k
        real(real64) :: initial_norm
        real(real64) :: threshold
        real(real64) :: trial_norm
        real(real64) :: alpha
        integer :: ndof
        integer :: nfree
        integer :: nc
        integer :: n
        integer :: i
        integer :: iteration
        integer :: backtrack
        integer :: lapack_info
        logical :: linear_ok
        logical :: accepted
        logical :: valid

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
        nc = mesh%number_constraints()
        if (nc > 0 .and. (opts%sparse .or. nfree == 0)) error stop DYN_INVALID_INPUT_ERROR
        allocate(lam(nc), source=0.0d0)

        ! Work on a private candidate so failure never publishes a partial step.
        q = displacement
        q(constrained) = conditions%value
        if (nfree == 0) then
            status%converged = .true.
            displacement = q
            if (present(multipliers)) multipliers = lam
            if (present(info)) info = status
            return
        end if
        call evaluate(q, lam, residual, valid)
        status%residual_norm = norm2(residual)
        if (.not.valid) then
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
                n = nfree + nc
                if (allocated(augmented)) deallocate(augmented)
                allocate(augmented(n,n), source=0.0d0)
                augmented(1:nfree,1:nfree) = dense_k(free,free)
                if (nc > 0) then
                    call mesh%assemble_constraints(q, gcon, jcon)
                    augmented(nfree+1:n,1:nfree) = jcon(:,free)
                    augmented(1:nfree,nfree+1:n) = transpose(jcon(:,free))
                end if
                linear_ok = all(ieee_is_finite(augmented))
                if (linear_ok) then
                    ! Right-hand side [R_u; -g] for the system [K G^T; G 0].
                    rhs = residual
                    rhs(nfree+1:n) = -residual(nfree+1:n)
                    ! LU overwrites the matrix; retain it to check the true residual.
                    original_k = augmented
                    allocate(pivots(n))
                    call DGETRF(n, n, augmented, n, pivots, lapack_info)
                    linear_ok = lapack_info == 0
                    if (linear_ok) then
                        delta = rhs
                        call DGETRS('N', n, 1, augmented, n, pivots, delta, n, lapack_info)
                        linear_ok = lapack_info == 0 .and. all(ieee_is_finite(delta))
                        if (linear_ok) linear_ok = norm2(rhs - matmul(original_k, delta)) <= &
                            1.0d-10 * norm2(rhs)
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
                trial(free) = q(free) + alpha * delta(1:nfree)
                lam_trial = lam + alpha * delta(nfree+1:)
                if (all(ieee_is_finite(trial)) .and. all(ieee_is_finite(lam_trial))) then
                    call evaluate(trial, lam_trial, trial_residual, valid)
                    if (valid) then
                        trial_norm = norm2(trial_residual)
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
            lam = lam_trial
            residual = trial_residual
            status%residual_norm = trial_norm
            status%iterations = iteration
            if (present(progress)) then
                report%iteration = iteration
                report%residual_norm = trial_norm
                report%tolerance = threshold
                report%step_fraction = alpha
                call progress(report, args = args)
            end if
        end do
        if (status%residual_norm > threshold) then
            call failure()
            return
        end if
        status%converged = .true.
        displacement = q
        if (present(multipliers)) multipliers = lam
        if (present(info)) info = status
    contains
        subroutine evaluate(qq, lamv, rv, is_valid)
            !! Evaluate the free residual `[f - f_int - G^T lambda; g]` and
            !! report whether all of its entries are finite.
            real(real64), intent(in), dimension(:) :: qq
                !! Full displacement vector.
            real(real64), intent(in), dimension(:) :: lamv
                !! Current Lagrange multipliers.
            real(real64), allocatable, intent(out), dimension(:) :: rv
                !! Free momentum residual followed by constraint values.
            logical, intent(out) :: is_valid
                !! True when every entry of `rv` is finite.

            real(real64), allocatable, dimension(:) :: fi
            real(real64), allocatable, dimension(:) :: gv
            real(real64), allocatable, dimension(:,:) :: jv

            call mesh%assemble_internal_force(qq, fi, nthreads=opts%nthreads)
            allocate(rv(nfree + nc), source=0.0d0)
            rv(1:nfree) = force(free) - fi(free)
            if (nc > 0) then
                call mesh%assemble_constraints(qq, gv, jv)
                rv(nfree+1:) = gv
                rv(1:nfree) = rv(1:nfree) - matmul(transpose(jv(:,free)), lamv)
            end if
            is_valid = all(ieee_is_finite(rv))
        end subroutine

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

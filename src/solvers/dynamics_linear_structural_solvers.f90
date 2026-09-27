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
module dynamics_linear_structural_solvers
    !! Advances a linear structural system with the generalized-alpha method.
    !! Dense and CSR integrators solve the semidiscrete equation
    !!
    !! $$M a(t) + C v(t) + K u(t) = f(t).$$
    !!
    !! Matrices and loads must already reflect any prescribed boundary
    !! conditions. The caller supplies a consistent initial acceleration and
    !! advances the load at both ends of each time step. Initialized objects
    !! reuse their work arrays and their effective-system factorization when
    !! the time step is unchanged. The linalg GMRES routine still allocates
    !! its own iteration workspace for large CSR systems.
    use iso_fortran_env, only : int32, real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use dynamics_error_handling
    use linalg, only : csr_matrix, msr_matrix, size, matmul, &
        csr_to_dense, pgmres_solver, operator(+), operator(*)
    use lapack, only : DGETRF, DGETRS
    implicit none
    private
    public :: structural_integrator
    public :: dense_generalized_alpha_integrator
    public :: sparse_generalized_alpha_integrator
    
    type, abstract :: structural_integrator
        !! Base class for linear structural time integrators.
        real(real64), private :: rho_infinity = 1.0d0
        real(real64), private :: alpha_m = 0.5d0
        real(real64), private :: alpha_f = 0.5d0
        real(real64), private :: gamma = 0.5d0
        real(real64), private :: beta = 0.25d0
        real(real64), private, allocatable, dimension(:) :: predicted_displacement
        real(real64), private, allocatable, dimension(:) :: predicted_velocity
        real(real64), private, allocatable, dimension(:) :: rhs
        real(real64), private, allocatable, dimension(:) :: next_acceleration
    contains
        procedure(integrator_step), deferred, public, pass :: step
            !! Advances the state by one time step.
        procedure, public :: solve => integrator_solve
            !! Advances a state through successive columns of a force history.
    end type

    interface
        subroutine integrator_step(this, force_current, force_next, dt, &
            displacement, velocity, acceleration)
            use iso_fortran_env, only : real64
            import structural_integrator
            class(structural_integrator), intent(inout) :: this
                !! The [[structural_integrator]] object.
            real(real64), intent(in), dimension(:) :: force_current
                !! External N-vector at the beginning of the step.
            real(real64), intent(in), dimension(:) :: force_next
                !! External N-vector at the end of the step.
            real(real64), intent(in) :: dt
                !! Positive step size \(h\).
            real(real64), intent(inout), dimension(:) :: displacement
                !! N-vector state at the start of the step, overwritten with 
                !! the end state.
            real(real64), intent(inout), dimension(:) :: velocity
                !! N-vector state at the start of the step, overwritten with 
                !! the end state.
            real(real64), intent(inout), dimension(:) :: acceleration
                !! N-vector state at the start of the step, overwritten with 
                !! the end state.
        end subroutine
    end interface

! ------------------------------------------------------------------------------
    type, extends(structural_integrator) :: dense_generalized_alpha_integrator
        !! Defines a generalized alpha integrator for dense matrix systems.
        real(real64), private, allocatable, dimension(:,:) :: mass
        real(real64), private, allocatable, dimension(:,:) :: damping
        real(real64), private, allocatable, dimension(:,:) :: stiffness
        real(real64), private, allocatable, dimension(:,:) :: lu
        integer(int32), private, allocatable, dimension(:) :: pivot
        real(real64), private :: cached_dt = -1.0d0
    contains
        procedure, public :: initialize => initialize_dense_integrator
            !! Copies system matrices and allocates reusable workspace.
        procedure, public :: step => dense_integrator_step
            !! Advances one step using cached dense LU when dt is unchanged.
    end type

! ------------------------------------------------------------------------------
    type, extends(structural_integrator) :: sparse_generalized_alpha_integrator
        !! Defines a generalized alpha integrator for sparse matrix systems.
        type(csr_matrix), private :: mass
        type(csr_matrix), private :: damping
        type(csr_matrix), private :: stiffness
        type(csr_matrix), private :: effective_matrix
        type(msr_matrix), private :: preconditioner
        real(real64), private, allocatable, dimension(:,:) :: small_lu
        real(real64), private, allocatable, dimension(:) :: diagonal
        integer(int32), private, allocatable, dimension(:) :: pivot
        integer(int32), private, allocatable, dimension(:) :: upper_start
        real(real64), private :: cached_dt = -1.0d0
    contains
        procedure, public :: initialize => initialize_sparse_integrator
            !! Copies CSR matrices and allocates reusable workspace.
        procedure, public :: step => sparse_integrator_step
            !! Advances one step using a cached effective CSR system.
    end type

contains
! ******************************************************************************
! STRUCTURAL_INTEGRATOR
! ------------------------------------------------------------------------------
subroutine integrator_solve(this, forces, dt, displacement, velocity, acceleration)
    !! Advances through a force history with one column per time point. The
    !! number of steps is size(forces, 2) - 1; only the final state is returned.
    class(structural_integrator), intent(inout) :: this
        !! Initialized dense or sparse integrator.
    real(real64), intent(in), dimension(:,:) :: forces
        !! N-by-(number of steps + 1) array of external forces.
    real(real64), intent(in) :: dt
        !! Positive constant time step.
    real(real64), intent(inout), dimension(:) :: displacement
        !! Initial state on input, final state on output.
    real(real64), intent(inout), dimension(:) :: velocity
        !! Initial state on input, final state on output.
    real(real64), intent(inout), dimension(:) :: acceleration
        !! Initial state on input, final state on output.
    integer(int32) :: step_index

    if (.not.allocated(this%rhs)) error stop DYN_INVALID_INPUT_ERROR
    if (size(forces, 2) < 2) error stop DYN_ARRAY_SIZE_ERROR
    if (size(forces, 1) /= size(displacement)) error stop DYN_ARRAY_SIZE_ERROR
    do step_index = 1, size(forces, 2) - 1
        call this%step(forces(:,step_index), forces(:,step_index+1), dt, &
            displacement, velocity, acceleration)
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine initialize_workspace(this, n, spectral_radius)
    !! Initializes the workspace for the integrator.
    class(structural_integrator), intent(inout) :: this
        !! The [[structural_integrator]] object.
    integer(int32), intent(in) :: n
        !! The number of degrees-of-freedom.
    real(real64), intent(in), optional :: spectral_radius
        !! High-frequency spectral radius in [0, 1], default 1.

    ! Local Variables
    real(real64) :: rho

    ! Process
    rho = 1.0d0
    if (present(spectral_radius)) rho = spectral_radius
    if (.not.ieee_is_finite(rho) .or. rho < 0.0d0 .or. rho > 1.0d0) &
        error stop DYN_INVALID_INPUT_ERROR

    this%rho_infinity = rho
    this%alpha_m = (2.0d0 * rho - 1.0d0) / (rho + 1.0d0)
    this%alpha_f = rho / (rho + 1.0d0)
    this%gamma = 0.5d0 + this%alpha_f - this%alpha_m
    this%beta = 0.25d0 * (1.0d0 + this%alpha_f - this%alpha_m)**2

    if (allocated(this%rhs)) then
        deallocate( &
            this%predicted_displacement, &
            this%predicted_velocity, &
            this%rhs, &
            this%next_acceleration &
        )
    end if
    allocate( &
        this%predicted_displacement(n), &
        this%predicted_velocity(n), &
        this%rhs(n), &
        this%next_acceleration(n) &
    )
end subroutine

! ******************************************************************************
! DENSE_GENERALIZED_ALPHA_INTEGRATOR
! ------------------------------------------------------------------------------
subroutine initialize_dense_integrator(this, mass, damping, stiffness, rho_infinity)
    !! Copies dense matrices and computes generalized-alpha coefficients once.
    class(dense_generalized_alpha_integrator), intent(inout) :: this
        !! Integrator to initialize or reinitialize.
    real(real64), intent(in), dimension(:,:) :: mass
        !! The N-by-N mass matrix.
    real(real64), intent(in), dimension(:,:) :: damping
        !! The N-by-N damping matrix.
    real(real64), intent(in), dimension(:,:) :: stiffness
        !! The N-by-N stiffness matrix.
    real(real64), intent(in), optional :: rho_infinity
        !! High-frequency spectral radius in [0, 1], default 1.

    ! Local Variables
    integer(int32) :: n

    ! Process
    n = size(mass, 1)
    if (n < 1 .or. size(mass, 2) /= n .or. &
        size(damping, 1) /= n .or. size(damping, 2) /= n .or. &
        size(stiffness, 1) /= n .or. size(stiffness, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR
    
    call initialize_workspace(this, n, rho_infinity)
    if (allocated(this%mass)) then
        deallocate( &
            this%mass, &
            this%damping, &
            this%stiffness, &
            this%lu, &
            this%pivot &
        )
    end if
    allocate(this%lu(n,n), this%pivot(n))
    allocate(this%mass(n, n), source = mass)
    allocate(this%damping(n, n), source = damping)
    allocate(this%stiffness(n, n), source = stiffness)
    this%cached_dt = -1.0d0
end subroutine

! ------------------------------------------------------------------------------
subroutine dense_integrator_step(this, force_current, force_next, dt, &
    displacement, velocity, acceleration)
    !! Advances a linear structural system by one generalized-alpha step.
    !! Let the input state be indexed by n and the output state by n+1. The
    !! method solves the weighted equilibrium
    !!
    !! $$M[(1-\alpha_m)a_{n+1}+\alpha_m a_n]
    !! + C[(1-\alpha_f)v_{n+1}+\alpha_f v_n]
    !! + K[(1-\alpha_f)u_{n+1}+\alpha_f u_n]
    !! = (1-\alpha_f)f_{n+1}+\alpha_f f_n.$$
    !!
    !! The high-frequency spectral radius \(\rho_\infty\) determines
    !!
    !! $$\alpha_m=\frac{2\rho_\infty-1}{\rho_\infty+1},\qquad
    !! \alpha_f=\frac{\rho_\infty}{\rho_\infty+1},\qquad
    !! \gamma=\tfrac12+\alpha_f-\alpha_m,\qquad
    !! \beta=\tfrac14(1+\alpha_f-\alpha_m)^2.$$
    !!
    !! With \(u^*=u_n+h v_n+h^2(\tfrac12-\beta)a_n\) and
    !! \(v^*=v_n+h(1-\gamma)a_n\), solve \(A a_{n+1}=b\), where
    !!
    !! $$A=(1-\alpha_m)M+(1-\alpha_f)\gamma h C
    !! +(1-\alpha_f)\beta h^2 K,$$
    !!
    !! $$b=(1-\alpha_f)f_{n+1}+\alpha_f f_n-\alpha_m M a_n
    !! -C[(1-\alpha_f)v^*+\alpha_f v_n]
    !! -K[(1-\alpha_f)u^*+\alpha_f u_n].$$
    !!
    !! Finally, \(u_{n+1}=u^*+\beta h^2 a_{n+1}\) and
    !! \(v_{n+1}=v^*+\gamma h a_{n+1}\). At \(\rho_\infty=1\), this is
    !! the average-acceleration (trapezoidal) Newmark scheme. 
    !! 
    !! The LU factorization of A is reused whenever the time step is unchanged.
    class(dense_generalized_alpha_integrator), intent(inout) :: this
        !! The [[dense_generalized_alpha_integrator]] object.
    real(real64), intent(in), dimension(:) :: force_current
        !! The current N-element external forcing vector.
    real(real64), intent(in), dimension(:) :: force_next
        !! The N-element external forcing vector at t + dt.
    real(real64), intent(in) :: dt
        !! The time step.
    real(real64), intent(inout), dimension(:) :: displacement
        !! The N-element displacement state vector.  On output, this vector is
        !! updated to the state at t + dt.
    real(real64), intent(inout), dimension(:) :: velocity
        !! The N-element velocity state vector.  On output, this vector is
        !! updated to the state at t + dt.
    real(real64), intent(inout), dimension(:) :: acceleration
        !! The N-element acceleration state vector.  On output, this vector is
        !! updated to the state at t + dt.

    ! Local Variables
    integer(int32) :: n, info

    ! Input Checking
    if (.not.allocated(this%mass)) error stop DYN_INVALID_INPUT_ERROR
    n = size(this%rhs)
    if (size(displacement) /= n .or. size(velocity) /= n .or. &
        size(acceleration) /= n .or. size(force_current) /= n .or. &
        size(force_next) /= n) error stop DYN_ARRAY_SIZE_ERROR
    if (.not.ieee_is_finite(dt) .or. dt <= 0.0d0) error stop DYN_INVALID_INPUT_ERROR

    ! A new dt changes the effective matrix; refactor only in that case.
    if (dt /= this%cached_dt) then
        this%lu = (1.0d0 - this%alpha_m) * this%mass + &
            (1.0d0 - this%alpha_f) * this%gamma * dt * this%damping + &
            (1.0d0 - this%alpha_f) * this%beta * dt**2 * this%stiffness
        call DGETRF(n, n, this%lu, n, this%pivot, info)
        if (info /= 0) error stop DYN_CONVERGENCE_ERROR
        this%cached_dt = dt
    end if

    ! Compute the predictor
    this%predicted_displacement = displacement + dt * velocity + &
        dt**2 * (0.5d0 - this%beta) * acceleration
    this%predicted_velocity = velocity + dt * (1.0d0 - this%gamma) * acceleration

    ! Compute the right-hand-side of the linear system
    this%rhs = (1.0d0 - this%alpha_f) * force_next + this%alpha_f * force_current - &
        this%alpha_m * matmul(this%mass, acceleration) - &
        matmul(this%damping, (1.0d0 - this%alpha_f) * this%predicted_velocity + &
            this%alpha_f * velocity) - &
        matmul(this%stiffness, (1.0d0 - this%alpha_f) * this%predicted_displacement + &
            this%alpha_f * displacement)
    this%next_acceleration = this%rhs

    ! Solve the linear system and update the state vectors
    call DGETRS('N', n, 1, this%lu, n, this%pivot, this%next_acceleration, n, info)
    if (info /= 0) error stop DYN_CONVERGENCE_ERROR
    displacement = this%predicted_displacement + this%beta * dt**2 * this%next_acceleration
    velocity = this%predicted_velocity + this%gamma * dt * this%next_acceleration
    acceleration = this%next_acceleration
end subroutine

! ******************************************************************************
! SPARSE_GENERALIZED_ALPHA_INTEGRATOR
! ------------------------------------------------------------------------------
subroutine initialize_sparse_integrator(this, mass, damping, stiffness, rho_infinity)
    !! Copies CSR matrices and prepares sparse or small-system dense workspace.
    class(sparse_generalized_alpha_integrator), intent(inout) :: this
        !! Integrator to initialize or reinitialize.
    type(csr_matrix), intent(in) :: mass
        !! The N-by-N mass matrix.
    type(csr_matrix), intent(in) :: damping
        !! The N-by-N damping matrix.
    type(csr_matrix), intent(in) :: stiffness
        !! The N-by-N stiffness matrix.
    real(real64), intent(in), optional :: rho_infinity
        !! High-frequency spectral radius in [0, 1], default 1.

    ! Local Variables
    integer(int32) :: n

    ! Process
    n = size(mass, 1)
    if (n < 1 .or. size(mass, 2) /= n .or. &
        size(damping, 1) /= n .or. size(damping, 2) /= n .or. &
        size(stiffness, 1) /= n .or. size(stiffness, 2) /= n) error stop DYN_MATRIX_SIZE_ERROR
    call initialize_workspace(this, n, rho_infinity)
    this%mass = mass
    this%damping = damping
    this%stiffness = stiffness
    if (allocated(this%small_lu)) deallocate(this%small_lu, this%pivot)
    if (allocated(this%diagonal)) then
        deallocate( &
            this%diagonal, &
            this%upper_start, &
            this%preconditioner%values, &
            this%preconditioner%indices &
        )
    end if
    if (n <= 32) then
        allocate(this%small_lu(n,n), this%pivot(n))
    else
        allocate( &
            this%diagonal(n), &
            this%upper_start(n), &
            this%preconditioner%values(n+1), &
            this%preconditioner%indices(n+1) &
        )
        this%preconditioner%m = n
        this%preconditioner%n = n
        this%preconditioner%nnz = 0
        this%preconditioner%values(n+1) = 0.0d0
        this%preconditioner%indices = n + 2
        this%upper_start = n + 2
    end if
    this%cached_dt = -1.0d0
end subroutine

! ------------------------------------------------------------------------------
subroutine sparse_integrator_step(this, force_current, force_next, dt, &
    displacement, velocity, acceleration)
    !! Advances by the generalized-alpha equations documented by
    !! dense_integrator_step using CSR matrix products. For N > 32,
    !! restarted GMRES uses a diagonal MSR preconditioner and checks the true
    !! residual; smaller systems use cached dense LU.
    class(sparse_generalized_alpha_integrator), intent(inout) :: this
        !! The [[sparse_generalized_alpha_integrator]] object.
    real(real64), intent(in), dimension(:) :: force_current
        !! The current N-element external forcing vector.
    real(real64), intent(in), dimension(:) :: force_next
        !! The N-element external forcing vector at t + dt.
    real(real64), intent(in) :: dt
        !! The time step.
    real(real64), intent(inout), dimension(:) :: displacement
        !! The N-element displacement state vector.  On output, this vector is
        !! updated to the state at t + dt.
    real(real64), intent(inout), dimension(:) :: velocity
        !! The N-element velocity state vector.  On output, this vector is
        !! updated to the state at t + dt.
    real(real64), intent(inout), dimension(:) :: acceleration
        !! The N-element acceleration state vector.  On output, this vector is
        !! updated to the state at t + dt.

    ! Local Variables
    integer(int32) :: n, info

    if (.not.allocated(this%rhs)) error stop DYN_INVALID_INPUT_ERROR
    n = size(this%rhs)
    if (size(displacement) /= n .or. size(velocity) /= n .or. &
        size(acceleration) /= n .or. size(force_current) /= n .or. &
        size(force_next) /= n) error stop DYN_ARRAY_SIZE_ERROR
    if (.not.ieee_is_finite(dt) .or. dt <= 0.0d0) error stop DYN_INVALID_INPUT_ERROR

    ! Refresh the effective system and its solver data only when dt changes.
    if (dt /= this%cached_dt) then
        this%effective_matrix = (1.0d0 - this%alpha_m) * this%mass + &
            (1.0d0 - this%alpha_f) * this%gamma * dt * this%damping + &
            (1.0d0 - this%alpha_f) * this%beta * dt**2 * this%stiffness
        if (n <= 32) then
            this%small_lu = csr_to_dense(this%effective_matrix)
            call DGETRF(n, n, this%small_lu, n, this%pivot, info)
            if (info /= 0) error stop DYN_CONVERGENCE_ERROR
        else
            this%diagonal = this%effective_matrix%extract_diagonal()
            if (any(.not.ieee_is_finite(this%diagonal)) .or. &
                any(abs(this%diagonal) <= tiny(1.0d0))) error stop DYN_CONVERGENCE_ERROR
            this%preconditioner%values(1:n) = 1.0d0 / this%diagonal
        end if
        this%cached_dt = dt
    end if

    this%predicted_displacement = displacement + dt * velocity + &
        dt**2 * (0.5d0 - this%beta) * acceleration
    this%predicted_velocity = velocity + dt * (1.0d0 - this%gamma) * acceleration
    this%rhs = (1.0d0 - this%alpha_f) * force_next + this%alpha_f * force_current - &
        this%alpha_m * matmul(this%mass, acceleration) - &
        matmul(this%damping, (1.0d0 - this%alpha_f) * this%predicted_velocity + &
            this%alpha_f * velocity) - &
        matmul(this%stiffness, (1.0d0 - this%alpha_f) * this%predicted_displacement + &
            this%alpha_f * displacement)
    if (n <= 32) then
        this%next_acceleration = this%rhs
        call DGETRS('N', n, 1, this%small_lu, n, this%pivot, this%next_acceleration, n, info)
        if (info /= 0) error stop DYN_CONVERGENCE_ERROR
    else
        this%next_acceleration = pgmres_solver(this%effective_matrix, &
            this%preconditioner, this%upper_start, this%rhs, &
            im = min(n, 100), tol = sqrt(epsilon(1.0d0)), maxits = 200)
        if (any(.not.ieee_is_finite(this%next_acceleration))) error stop DYN_CONVERGENCE_ERROR
        if (norm2(this%rhs - matmul(this%effective_matrix, this%next_acceleration)) > &
            1.0d-8 * max(norm2(this%rhs), 1.0d0)) error stop DYN_CONVERGENCE_ERROR
    end if
    displacement = this%predicted_displacement + this%beta * dt**2 * this%next_acceleration
    velocity = this%predicted_velocity + this%gamma * dt * this%next_acceleration
    acceleration = this%next_acceleration
end subroutine

! ------------------------------------------------------------------------------
end module
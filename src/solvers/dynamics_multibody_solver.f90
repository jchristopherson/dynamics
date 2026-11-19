! SPDX-License-Identifier: MIT
module dynamics_multibody_solver
    !! Fixed-step, dense generalized-alpha integration of spatial multibody
    !! models. Holonomic constraints are enforced at the new position and
    !! velocities are projected onto the endpoint tangent space.
    use iso_fortran_env, only : real64, int32
    use ieee_arithmetic, only : ieee_is_finite
    use linalg, only : lu_factor, solve_lu
    use dynamics_multibody, only : multibody_model
    implicit none
    private
    public :: multibody_integrator, multibody_state

    type :: multibody_state
        real(real64) :: time = 0.0d0
            !! Time at the position and velocity sample.
        real(real64), allocatable :: position(:), velocity(:)
            !! Assembled coordinates and their time derivatives.
        real(real64), allocatable :: acceleration(:)
            !! Newmark acceleration; after startup this is an integration variable.
        real(real64), allocatable :: multipliers(:)
            !! Startup or force-stage multipliers, with sign M*a + G^T*lambda = f.
    end type

    type :: multibody_integrator
        real(real64) :: rho_infinity = 0.8d0
            !! High-frequency spectral radius of the underlying generalized-alpha method, in [0,1].
        real(real64) :: force_tolerance = 1.0d-9
            !! Mixed absolute/relative Newton tolerance for the equations of motion.
        real(real64) :: constraint_tolerance = 1.0d-10
            !! Absolute infinity-norm tolerance for position constraint residuals.
        integer :: max_iterations = 30
            !! Maximum Newton corrections per fixed step.
    contains
        procedure :: initialize
        procedure :: step
        procedure :: solve
    end type

contains

subroutine initialize(this, model, state, time, position, velocity)
    !! Initializes from a position on the constraint manifold. The supplied
    !! velocity is mass-projected before solving for consistent acceleration.
    class(multibody_integrator), intent(in) :: this
    type(multibody_model), intent(in) :: model
    type(multibody_state), intent(out) :: state
    real(real64), intent(in) :: time, position(:), velocity(:)
    integer :: n, nc
    real(real64), allocatable :: c(:), g(:,:), m(:,:), force(:), answer(:)

    call check_settings(this)
    n = model%get_dof_count()
    nc = model%get_constraint_count()
    if (n == 0 .or. nc > n .or. size(position) /= n .or. size(velocity) /= n) &
        error stop "Invalid multibody initial state dimensions"
    if (.not. ieee_is_finite(time) .or. .not. all(ieee_is_finite(position)) .or. &
        .not. all(ieee_is_finite(velocity))) error stop "Nonfinite multibody initial state"
    c = model%constraints(time, position)
    if (.not. all(ieee_is_finite(c)) .or. infinity_norm(c) > this%constraint_tolerance) &
        error stop "Initial multibody position does not satisfy the constraints"

    state%time = time
    state%position = position
    state%velocity = velocity
    g = constraint_jacobian(model, time, position)
    m = model%mass_matrix(position)
    call project_velocity(model, time, position, m, g, state%velocity)
    allocate(force(n))
    call model%equations(time, position, state%velocity, m, force)
    answer = saddle_solve(m, g, force, -constraint_curvature(model, time, position, state%velocity))
    state%acceleration = answer(:n)
    state%multipliers = answer(n+1:)
end subroutine

subroutine step(this, model, state, dt)
    !! Advances an initialized state by a positive time step. A dense Newton
    !! solve enforces endpoint position, followed by a velocity projection.
    !! Invalid, rank-deficient or nonconvergent steps stop with a diagnostic.
    class(multibody_integrator), intent(in) :: this
    type(multibody_model), intent(in) :: model
    type(multibody_state), intent(inout) :: state
    real(real64), intent(in) :: dt
    real(real64) :: am, af, beta, gamma, bdt2, gdt, endpoint, stage_time, h, fraction, score, trial_score
    real(real64) :: force_limit, position_limit
    real(real64), allocatable :: x(:), trial(:), residual(:), trial_residual(:), correction(:), jac(:,:)
    real(real64), allocatable :: qpred(:), vpred(:), qnew(:), g(:,:), stage_g(:,:), m(:,:)
    integer :: n, nc, iteration, j, search

    call check_settings(this)
    n = model%get_dof_count()
    nc = model%get_constraint_count()
    if (.not. allocated(state%position) .or. .not. allocated(state%velocity) .or. &
        .not. allocated(state%acceleration) .or. .not. allocated(state%multipliers)) &
        error stop "Initialize the multibody state before stepping"
    if (n == 0 .or. nc > n .or. size(state%position) /= n .or. size(state%velocity) /= n .or. &
        size(state%acceleration) /= n .or. size(state%multipliers) /= nc) error stop "Invalid multibody state dimensions"
    if (.not. ieee_is_finite(dt) .or. dt <= 0.0d0) error stop "Invalid multibody step size"
    endpoint = state%time + dt
    if (.not. ieee_is_finite(endpoint) .or. endpoint <= state%time) error stop "Unrepresentable multibody step time"

    am = (2*this%rho_infinity-1)/(1+this%rho_infinity)
    af = this%rho_infinity/(1+this%rho_infinity)
    gamma = 0.5d0 + af - am
    beta = 0.25d0*(1+af-am)**2
    bdt2 = beta*dt**2
    gdt = gamma*dt
    if (.not. ieee_is_finite(bdt2) .or. bdt2 <= tiny(bdt2)) error stop "Unrepresentable multibody step size"
    stage_time = state%time + (1-af)*dt
    qpred = state%position + dt*state%velocity + (0.5d0-beta)*dt**2*state%acceleration
    vpred = state%velocity + (1-gamma)*dt*state%acceleration
    x = [state%acceleration, state%multipliers]
    allocate(jac(n+nc,n+nc), stage_g(nc,n))
    residual = evaluate(x, stage_g)
    force_limit = this%force_tolerance*max(1.0d0,infinity_norm(residual(:n)))
    ! Minimize constraint drift without asking for sub-roundoff coordinate changes.
    position_limit = min(this%constraint_tolerance, max(this%force_tolerance*bdt2, &
        16*epsilon(1.0d0)*max(1.0d0,infinity_norm(qpred))))
    score = merit(residual)
    do iteration = 1, this%max_iterations
        if (score <= 1.0d0) exit
        jac = 0.0d0
        do j = 1, n
            h = epsilon(1.0d0)**0.25d0 * max(1.0d0,abs(x(j)))
            trial = x
            trial(j) = x(j) + h
            trial_residual = evaluate(trial)
            trial(j) = x(j) - h
            correction = evaluate(trial)
            jac(:n,j) = (trial_residual(:n)-correction(:n))/(2*h)
        end do
        qnew = qpred + bdt2*x(:n)
        jac(n+1:,:n) = constraint_jacobian(model, endpoint, qnew)
        jac(:n,n+1:) = transpose(stage_g)
        correction = linear_solve(jac, -residual)
        fraction = 1.0d0
        do search = 1, 12
            trial = x + fraction*correction
            trial_residual = evaluate(trial, stage_g)
            trial_score = merit(trial_residual)
            if (trial_score < score) exit
            fraction = 0.5d0*fraction
        end do
        if (trial_score >= score) error stop "Multibody Newton line search stalled; reduce the time step"
        x = trial
        residual = trial_residual
        score = trial_score
    end do
    if (score > 1.0d0) error stop "Multibody Newton iteration did not converge"

    qnew = qpred + bdt2*x(:n)
    g = constraint_jacobian(model, endpoint, qnew)
    m = model%mass_matrix(qnew)
    state%velocity = vpred + gdt*x(:n)
    call project_velocity(model, endpoint, qnew, m, g, state%velocity)
    state%time = endpoint
    state%position = qnew
    state%acceleration = x(:n)
    state%multipliers = x(n+1:)
contains
    function evaluate(y, gradient) result(r)
        !! Unknowns are the endpoint acceleration and the force-stage
        !! multipliers. Divide position rows by beta*dt**2 to balance Newton.
        real(real64), intent(in) :: y(:)
        real(real64), intent(out), optional :: gradient(:,:)
        real(real64) :: r(n+nc), q(n), v(n), a(n), mass(n,n), forces(n), gs(nc,n), qe(n)

        qe = qpred + bdt2*y(:n)
        q = (1-af)*qe + af*state%position
        v = (1-af)*(vpred + gdt*y(:n)) + af*state%velocity
        a = (1-am)*y(:n) + am*state%acceleration
        call model%equations(stage_time, q, v, mass, forces)
        gs = constraint_jacobian(model, stage_time, q)
        if (present(gradient)) gradient = gs
        r(:n) = matmul(mass,a) + matmul(transpose(gs),y(n+1:)) - forces
        r(n+1:) = model%constraints(endpoint, qe)/bdt2
    end function

    real(real64) function merit(r) result(value)
        !! Line-search norm tests force and position convergence independently.
        real(real64), intent(in) :: r(:)
        value = huge(1.0d0)
        if (.not. all(ieee_is_finite(r))) return
        value = max(infinity_norm(r(:n))/force_limit, infinity_norm(r(n+1:))*bdt2/position_limit)
    end function
end subroutine

function solve(this, model, state, dt, nstep) result(history)
    !! Returns nstep+1 samples including the initialized state; updates state
    !! in place to the last sample. Use step directly to avoid storing history.
    class(multibody_integrator), intent(in) :: this
    type(multibody_model), intent(in) :: model
    type(multibody_state), intent(inout) :: state
    real(real64), intent(in) :: dt
    integer, intent(in) :: nstep
    type(multibody_state), allocatable :: history(:)
    integer :: i

    if (nstep < 0) error stop "Invalid multibody step count"
    allocate(history(nstep+1))
    history(1) = state
    do i = 1, nstep
        call this%step(model, state, dt)
        history(i+1) = state
    end do
end function

function constraint_jacobian(model, time, q) result(g)
    !! Centered coordinate derivative of the assembled holonomic residuals.
    type(multibody_model), intent(in) :: model
    real(real64), intent(in) :: time, q(:)
    real(real64) :: g(model%get_constraint_count(),size(q)), qp(size(q)), qm(size(q)), h
    integer :: j

    if (size(g,1) == 0) return
    do j = 1, size(q)
        h = epsilon(1.0d0)**(1.0d0/3)*max(1.0d0,abs(q(j)))
        qp = q
        qm = q
        qp(j) = qp(j) + h
        qm(j) = qm(j) - h
        g(:,j) = (model%constraints(time,qp)-model%constraints(time,qm))/(2*h)
    end do
end function

function constraint_curvature(model, time, q, v) result(c)
    !! Second derivative of Phi(t+s,q+s*v) at s=0; the missing G*a
    !! term is obtained from the startup saddle-point solve.
    type(multibody_model), intent(in) :: model
    real(real64), intent(in) :: time, q(:), v(:)
    real(real64) :: c(model%get_constraint_count()), scale, h

    if (size(c) == 0) return
    scale = max(1.0d0,infinity_norm(v))
    h = epsilon(1.0d0)**0.25d0/scale
    c = (model%constraints(time+h,q+h*v) + model%constraints(time-h,q-h*v) - &
        2*model%constraints(time,q))/h**2
end function

subroutine project_velocity(model, time, q, m, g, v)
    !! Minimum kinetic-energy correction enforcing G*v + Phi_t = 0.
    type(multibody_model), intent(in) :: model
    real(real64), intent(in) :: time, q(:), m(:,:), g(:,:)
    real(real64), intent(inout) :: v(:)
    real(real64) :: h, defect(size(g,1)), zero(size(v)), answer(size(v)+size(g,1))

    if (size(defect) == 0) return
    h = epsilon(1.0d0)**(1.0d0/3)*max(1.0d0,abs(time))
    defect = matmul(g,v) + (model%constraints(time+h,q)-model%constraints(time-h,q))/(2*h)
    zero = 0.0d0
    answer = saddle_solve(m, g, zero, -defect)
    v = v + answer(:size(v))
end subroutine

function saddle_solve(m, g, force, constraint) result(answer)
    !! Solves the mechanical saddle system with a zero multiplier block.
    real(real64), intent(in) :: m(:,:), g(:,:), force(:), constraint(:)
    real(real64) :: answer(size(force)+size(constraint)), block(size(answer),size(answer))
    integer :: n

    n = size(force)
    block = 0.0d0
    block(:n,:n) = m
    block(:n,n+1:) = transpose(g)
    block(n+1:,:n) = g
    answer = linear_solve(block, [force,constraint])
end function

function linear_solve(a, rhs) result(x)
    !! Equilibrates the dense system before pivoted LU and diagnoses singular
    !! models (e.g., redundant constraints or unrestrained massless modes).
    real(real64), intent(in) :: a(:,:), rhs(:)
    real(real64) :: x(size(rhs)), work(size(rhs),size(rhs)), row(size(rhs)), col(size(rhs))
    real(real64), allocatable :: lu(:,:)
    integer(int32), allocatable :: pivot(:)
    integer :: i

    if (.not. all(ieee_is_finite(a)) .or. .not. all(ieee_is_finite(rhs))) error stop "Nonfinite multibody equations"
    row = maxval(abs(a),dim=2)
    if (any(row <= tiny(1.0d0))) error stop "Singular multibody system (zero row)"
    work = a / spread(row,2,size(rhs))
    col = maxval(abs(work),dim=1)
    if (any(col <= tiny(1.0d0))) error stop "Singular multibody system (zero column)"
    work = work / spread(col,1,size(rhs))
    call lu_factor(work, ipvt=pivot, lu=lu)
    do i = 1, size(rhs)
        if (.not. (abs(lu(i,i)) > 100*epsilon(1.0d0))) error stop "Rank-deficient multibody system"
    end do
    x = solve_lu(lu,pivot,rhs/row)/col
    if (.not. all(ieee_is_finite(x))) error stop "Nonfinite multibody system solution"
end function

pure real(real64) function infinity_norm(x) result(value)
    !! The infinity norm is zero for an empty (unconstrained) residual vector.
    real(real64), intent(in) :: x(:)
    value = max(0.0d0,maxval(abs(x)))
end function

subroutine check_settings(this)
    !! Rejects invalid controls before entering the numerical routines.
    class(multibody_integrator), intent(in) :: this
    if (.not. ieee_is_finite(this%rho_infinity) .or. this%rho_infinity < 0 .or. this%rho_infinity > 1) &
        error stop "Invalid multibody rho_infinity"
    if (.not. ieee_is_finite(this%force_tolerance) .or. this%force_tolerance <= 0 .or. &
        .not. ieee_is_finite(this%constraint_tolerance) .or. this%constraint_tolerance <= 0 .or. &
        this%max_iterations < 1) error stop "Invalid multibody Newton settings"
end subroutine

end module

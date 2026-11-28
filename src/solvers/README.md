## Solver Modules

The solver functionality in this directory covers structural equilibrium and time integration:

* **`dynamics_linear_structural_solvers`** advances the linear system \(M a + C v + K u = f\) with the generalized-alpha method. It provides dense and CSR sparse integrators, accepts loads at both ends of a time step, and reuses the effective-system factorization when the time step is unchanged. Apply prescribed boundary conditions before integration and provide a consistent initial acceleration.
* **`dynamics_variational_integrators`** advances constrained rigid-body systems in maximal coordinates using discrete Euler-Lagrange equations and Lagrange multipliers. It supports dense and graph-factorized Newton solvers, user force and holonomic-constraint callbacks, optional constraint Jacobians, and configurable force sampling at the left endpoint, implicit endpoint, or midpoint.
* **`dynamics_flexible_body_solvers`** solves nonlinear static equilibrium and transient dynamics of flexible meshes, with dense methods and a sparse BDF option for unconstrained and joint-constrained motion.

For transient setup, solver restrictions, and sparse scaling limitations, see
the [flexible-body transient guide](../flexible_bodies/README.md#transient-analysis).

## Nonlinear flexible static equilibrium

`solve_flexible_static` solves `internal_force(q) = force` at free DOFs.
Supply a full force vector and initial displacement in mesh equation order,
and `nodal_condition` records for prescribed displacements. Prescribed values
are imposed before evaluating the residual; external loads at constrained DOFs
do not affect equilibrium. To obtain reactions, assemble the internal force
after convergence and subtract the applied force at constrained DOFs.

```fortran
type(flexible_static_settings) :: settings
type(flexible_static_info) :: info

settings%sparse = .true.
call solve_flexible_static(mesh, force, supports, q, settings, info)
if (.not.info%converged) error stop "Static equilibrium failed"
```

The solver uses the element tangent, full Newton corrections, and residual-norm
backtracking. Convergence requires the free residual norm to be no greater than
`absolute_tolerance + relative_tolerance * max(initial_residual_norm, norm2(force_free))`.
The norm combines all free components without per-DOF scaling; models with
different force and moment scales may need nondimensionalization. Set tolerances
in units appropriate to the model. Dense tangents use LU;
CSR tangents use ILUT-preconditioned restarted GMRES with a checked true linear
residual. ILUT uses a drop tolerance of the square root of machine epsilon and
is rebuilt for each Newton tangent. The sparse method may still fail on poorly
conditioned or indefinite systems even when a solution exists.
`nthreads` has the same element thread-safety requirements as mesh assembly.
With `info`, convergence or linear-solver failure leaves `q` unchanged;
without it, failure triggers `error stop`. Invalid inputs also stop. With no
free DOFs, imposing the prescribed displacements succeeds without iteration.

Loads are fixed global forces, not follower loads; the solver does not implement
load stepping, arc-length continuation, contact, or time integration. Repeated
calls can apply successive loads or prescribed displacements, using the previous
solution as the next initial guess. Force and tangent evaluations must be
consistent and safe to repeat at trial displacements; the existing material
contract provides energy-based nonlinear laws, not commit/rollback for
history-dependent plasticity. Geometric and constitutive fidelity remain those
of the selected element and material.

## Progress reporting

Each solver accepts an optional progress callback, so reporting is off unless
you supply one. Callbacks are declared in
[`dynamics_solver_reporting`](./dynamics_solver_reporting.f90), except the
variational integrator callback, which is declared in its own module because it
receives a `variational_state`. Each callback receives an optional `args`
argument, which the solver forwards unchanged.

* **`solve_flexible_static`** reports each accepted Newton iterate through a
  `static_solver_callback`. The `static_solver_progress` argument contains the
  iteration count, free residual norm, convergence tolerance, and the fraction of
  the Newton correction accepted after backtracking. Rejected trial steps are
  not reported.
* **`integrator_solve`** (the `solve` method of the structural integrators)
  reports each completed time step through a `time_history_solver_callback`. The
  `time_history_solver_progress` argument contains the step number, total steps,
  elapsed time, and copies of the displacement, velocity, and acceleration.
* **`vi_solve`** reports each completed time step through a
  `variational_step_callback` that receives the step number, total steps, and the
  new `variational_state`. It does not report iteration counts or look-ahead
  steps.

Pass a procedure whose interface matches the callback:

```fortran
call integrator%solve(force_history, dt, displacement, velocity, acceleration, print_step)
```

## Variational Integrator Example
The [`variational_integrator_example`](../../examples/variational_integrator_example.f90) simulates a planar double pendulum in maximal coordinates. Both connecting rods have distributed mass, finite cross-section inertia, and gravity loading at their centers of mass. Six holonomic constraints pin the first rod to ground and join the two rod endpoints.

The example selects the graph-factorized solver from Brüdigam et al. (2023), supplies force and constraint callbacks, and provides an analytic reduced constraint Jacobian for efficient Newton iterations:

```fortran
type(rigid_body), dimension(2) :: bodies
type(variational_state) :: initial_state
type(variational_state), allocatable, dimension(:) :: solution
type(variational_integrator) :: integrator

! Each connecting rod carries its own mass and center-of-mass inertia tensor.
bodies(1) = rigid_body(mass1, rod_inertia(mass1, length1, width1))
bodies(2) = rigid_body(mass2, rod_inertia(mass2, length2, width2))

call initialize_variational_state(initial_state, 2)
initial_state%orientation(1) = quaternion(angle1_initial, &
    [0.0d0, 0.0d0, 1.0d0])
initial_state%orientation(2) = quaternion(angle2_initial, &
    [0.0d0, 0.0d0, 1.0d0])

! Set compatible center-of-mass positions for the two endpoint constraints.
direction1 = [sin(angle1_initial), -cos(angle1_initial), 0.0d0]
direction2 = [sin(angle2_initial), -cos(angle2_initial), 0.0d0]
initial_state%position(:,1) = 0.5d0 * length1 * direction1
initial_state%position(:,2) = length1 * direction1 + &
    0.5d0 * length2 * direction2

integrator%settings%linear_solver = VI_GRAPH_FACTORIZED_SOLVER
solution = integrator%solve(bodies, initial_state, dt, ntime, &
    constraint_count = 6, &
    constraint = pendulum_constraints, &
    force_function = gravity_forces, &
    constraint_jacobian = pendulum_constraint_jacobian, &
    args = parameters)
```

The state uses a position and orientation for each body's center of mass, so the endpoint constraints must be satisfied by the initial positions and orientations before solving. `constraint_count = 6` declares the scalar equations returned by the callback; the analytic Jacobian callback supplies their derivatives to Newton's method. The graph-factorized option exploits the sparse body/constraint connectivity of the system, while the dense solver is useful for small reference problems.

![Double-pendulum rod angles produced by the variational integrator example](../../images/variational_integrator_example.png?raw=true)
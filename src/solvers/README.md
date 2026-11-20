## Solver Modules

The solver functionality in this directory covers linear structural time integration and constrained rigid-body variational integration:

* **`dynamics_linear_structural_solvers`** advances the linear system \(M a + C v + K u = f\) with the generalized-alpha method. It provides dense and CSR sparse integrators, accepts loads at both ends of a time step, and reuses the effective-system factorization when the time step is unchanged. Apply prescribed boundary conditions before integration and provide a consistent initial acceleration.
* **`dynamics_variational_integrators`** advances constrained rigid-body systems in maximal coordinates using discrete Euler-Lagrange equations and Lagrange multipliers. It supports dense and graph-factorized Newton solvers, user force and holonomic-constraint callbacks, optional constraint Jacobians, and configurable force sampling at the left endpoint, implicit endpoint, or midpoint.

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
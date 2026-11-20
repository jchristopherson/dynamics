## Vibration Analysis Modules

* [`dynamics_vibrations`](dynamics_vibrations.f90) provides vibration and transient-response metrics, including quality factor, bandwidth, damping from logarithmic decrement or overshoot, and helpers for analyzing free and step responses.
* [`dynamics_stability`](dynamics_stability.f90) classifies fixed points from the eigenvalues of the linearized dynamics matrix. Eigenvalues with zero real part may require further nonlinear analysis.
* [`dynamics_maps`](dynamics_maps.f90) extracts Poincaré sections from either a sampled trajectory or an ODE model. It interpolates crossings between samples and can retain both crossing directions or only crossings from a specified side of the plane.

A Poincaré section records where a trajectory crosses a surface instead of plotting the entire time history. The examples below show how to use a geometric plane to slice an existing trajectory and how to sample a periodically forced oscillator at a fixed forcing phase.

## Lorenz System: Section of a Sampled Trajectory

[`poincare_example`](../../examples/poincare_example.f90) integrates the Lorenz equations and slices the trajectory at $z=25$. The `plane` constructor takes a point on the plane and a normal. The first column of `sol` holds time; the next three hold $x$, $y$, and $z$.

```fortran
pln = plane([0.0d0, 0.0d0, 25.0d0], [0.0d0, 0.0d0, 1.0d0])
mdl%fcn => eom
call integrator%solve(mdl, [0.0d0, 1.0d2], [1.0d0, 0.0d0, 1.0d0])
sol = integrator%get_solution()
map = poincare_map(sol(:,2), sol(:,3), sol(:,4), pln)
```

The call above uses the default two-sided section. Set `side = POINCARE_ONE_SIDED_FROM_FRONT` to keep crossings from $z>25$, or `side = POINCARE_ONE_SIDED_FROM_BACK` to keep crossings from $z<25$. The returned array has one row per retained intersection and three columns for its coordinates. Plot the first two columns for a view of the section, or overlay all three on the original trajectory.

| Crossing selection | Section in the $xy$ plane | Section overlaid on the trajectory |
| --- | --- | --- |
| Two-sided (default) | ![Lorenz section, both crossing directions](../../images/poincare_two_sided_2D.png) | ![Lorenz trajectory with two-sided section points](../../images/poincare_two_sided_3D.png) |
| From the front | ![Lorenz section, crossings from above](../../images/poincare_one_sided_front_2D.png) | ![Lorenz trajectory with front-side section points](../../images/poincare_one_sided_front_3D.png) |
| From the back | ![Lorenz section, crossings from below](../../images/poincare_one_sided_back_2D.png) | ![Lorenz trajectory with back-side section points](../../images/poincare_one_sided_back_3D.png) |

## Duffing Oscillator: Fixed-Phase Section

[`duffing_poincare_example`](../../examples/duffing_poincare_example.f90) models a damped, periodically forced Duffing oscillator. It adds forcing phase $\phi$ to the displacement and velocity state, with $\dot{\phi}=\omega$, so the forcing can be evaluated as $\gamma\cos\phi$. The section-coordinate callback returns `[state(1), state(2), sin(state(3))]`, so the default section plane selects zeros of $\sin\phi$. Keeping only crossings from the back selects phase zero modulo $2\pi$, rather than also sampling a half-period later.

```fortran
model%fcn => duffing_eom
call integrator%solve(model, [0.0d0, 0.5d0 * transient_end, &
    transient_end], [0.0d0, 0.0d0, 0.0d0])
transient_solution = integrator%get_solution()

section = poincare_map(model, [transient_end, response_end], &
    transient_solution(3,2:4), response_points, &
    side = POINCARE_ONE_SIDED_FROM_BACK, solver = integrator, &
    chunk_size = 20 * points_per_cycle, &
    coordinates = duffing_section_coordinates)
```

The example discards the first 50 forcing cycles, then analyzes 50,000 cycles at 50 sample intervals per cycle. The ODE overload processes the trajectory in chunks rather than storing the full history. The plot shows displacement versus velocity at the retained crossings; the markers are fixed-phase samples, not a continuous trajectory.

![Duffing oscillator displacement and velocity sampled once per forcing period](../../images/duffing_poincare_example.png)

To build the examples, configure CMake with `BUILD_DYNAMICS_EXAMPLES=ON` and build the `poincare_example` and `duffing_poincare_example` targets.
# DYNAMICS C API

The DYNAMICS C API is a C-compatible interface to the library. The public declarations are in [`dynamics.h`](../../src/c/dynamics.h).

## Conventions

- Matrices use column-major storage, matching the Fortran API. For a matrix with `m` rows, `n` columns, and leading dimension `ld`, element `(i, j)` is stored at `j * ld + i` using zero-based C indices.
- Mechanism link, joint, and frame indices are one-based, matching the Fortran API.
- Memory returned by an allocation or creation routine is owned by the caller and must be released with its corresponding `c_free_*` routine.
- An optional pointer may be passed as `NULL` only where the function declaration or reference explicitly permits it.
- Solver callbacks use C calling conventions and receive array dimensions as explicit arguments.
- Variational state histories use column-major `3 × nbody × ntime` arrays;
	multiplier histories use `nconstraint × ntime` arrays.
- The prescribed-motion linkage callback is synchronous and not reentrant.

## Linking

Include the installed header and link against the DYNAMICS library and its dependencies. A CMake consumer can use the exported DYNAMICS target after installation:

```cmake
find_package(dynamics CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE dynamics::dynamics)
```

The C interface is enabled when building DYNAMICS with:

```sh
cmake -S . -B build -DBUILD_DYNAMICS_C_INTERFACE=ON
```

## Example

The repository contains a complete closed-loop example at
[`examples/c_four_bar_example_1.c`](../../examples/c_four_bar_example_1.c). It shows
how to allocate mechanism links, connect them with `c_joint` values, create a
planar linkage, solve its forward kinematics, and release all allocated data.

## Beam coordinates

The beam element APIs use the local element coordinate system shown below. The
local x-axis runs from node 1 to node 2. For a 3D beam, the orientation point
defines the local z-axis and the local y-axis completes the right-handed frame.

![Beam element coordinate systems](beam_coordinate_system.svg)

## Denavit-Hartenberg parameters

The `c_dh_parameter_set` structure follows the standard Denavit-Hartenberg
convention used by the kinematics API. The transform maps coordinates from
frame `i` into frame `i-1` as

`Rz(joint_angle) Tz(link_offset) Tx(link_length) Rx(link_twist)`.

![Standard Denavit-Hartenberg parameters](dh_parameter_set.svg)

## Geometry operations

The geometry API provides point, line, plane, and Plucker-line representations,
along with projection, distance, parallelism, intersection, and common-normal
operations. The relationships are summarized below.

![Geometry representations and operations](geometry_operations.svg)

## Variational and linkage dynamics

The variational API supports two workflows:

1. `c_variational_integrator_solve` integrates caller-defined rigid bodies in
	maximal coordinates. C callbacks supply applied loads, equality constraints,
	and optionally an analytic reduced constraint Jacobian.
2. `c_create_serial_linkage_dynamic_model` and
	`c_create_linkage_dynamic_model` convert linkage descriptions into opaque
	dynamic-model handles. `c_linkage_dynamic_solve_joint` then applies gravity,
	constant body loads, and optional relative revolute or prismatic motion.
	`c_linkage_dynamic_solve` remains available for absolute planar body-angle input.

Initialize `c_variational_integrator_settings` with
`c_default_variational_integrator_settings` before changing individual fields.
The `force_evaluation` field selects the applied-load sampling rule for both
direct and linkage solves:

- `DYN_VI_FORCE_LEFT_ENDPOINT` (default) evaluates loads at the current state;
	this is explicit and least expensive.
- `DYN_VI_FORCE_IMPLICIT_ENDPOINT` evaluates loads at the trial next state in
	every Newton residual evaluation. This is more robust for stiff damping but
	increases nonlinear work and also endpoint-samples springs and user loads.
- `DYN_VI_FORCE_MIDPOINT` evaluates loads at the time/state midpoint in every
	Newton residual evaluation. Translation, time, and linear velocity are
	averaged; orientation uses the unit-quaternion geodesic midpoint, and angular
	velocity is averaged in world coordinates then expressed in the midpoint
	body frame. This is midpoint force sampling, not a complete implicit-midpoint
	state integrator.

For a scalar mass-damper equation `m * dv/dt = -c * v`, with $r=h c/m$, the velocity
amplification factors are $1-r$ (left endpoint), $1/(1+r)$ (implicit endpoint),
and $(1-r/2)/(1+r/2)$ (midpoint). Thus the explicit mode can grow energy for
$r>2$, the implicit mode strongly damps stiff modes, and midpoint is stable but
approaches a sign-alternating factor of $-1$ for very stiff modes. All applied
loads use the selected sampling rule, not only damper elements. Force callbacks
in the implicit and midpoint modes are called repeatedly during Newton and
should be deterministic functions of their input state, time, and user data.

The direct and linkage solve routines write caller-owned histories using these
column-major shapes:

- position, velocity, and angular velocity: `3 × nbody × ntime`;
- orientation: `nbody × ntime`;
- multipliers: `nconstraint × ntime`.

Multiplier column `i` belongs to simulation point `i`. The final column is
computed with a noncommitting look-ahead step. For prescribed linkage motion,
the driving multiplier is appended after all joint and planar constraints and
represents the required actuator torque (angular input) or axial force
(prismatic input), positive toward increasing prescribed coordinate.

Use `c_linkage_dynamic_joint_reactions` to convert one state and multiplier
column into world-frame joint forces and moments. Reactions are returned in
mechanism joint order and act on each joint's child link. The reaction on the
parent link is equal and opposite.

Linkage dynamic models may also contain linear force elements:

- `c_linear_spring` acts in tension and compression between arbitrary
	body-fixed points. Its `free_length` defines zero force and therefore any
	preload at the initial configuration.
- `c_linear_damper` acts only on relative velocity along the current element
	axis.
- `c_torsional_spring` acts about the axis of its referenced revolute joint;
	`free_angle` defines its zero-torque angle.
- `c_torsional_damper` opposes only relative twist rate about that joint axis.

Body index zero denotes a ground attachment whose point is expressed in world
coordinates. Other attachment points are expressed in their body frame. Use
the element count and result routines to query current length/rate/force or
angle/rate/torque values.

Opaque dynamic-model handles own copied linkage and mass-property data; release
them with `c_free_linkage_dynamic_model`. Callback state views are temporary and
must not be retained after a callback returns. The prescribed-motion callback
bridge is synchronous and not reentrant.

### Prescribed joint motion

`c_linkage_dynamic_solve_joint` has the same argument layout as the existing
`c_linkage_dynamic_solve`, but replaces `prescribed_body` with `prescribed_joint`.
The original function and its ABI are unchanged. Both use the scalar
`c_linkage_prescribed_motion` callback:

| Entry point / selected joint | Callback return value |
| --- | --- |
| `c_linkage_dynamic_solve_joint`, `DYN_REVOLUTE_JOINT` | Relative child-to-parent joint-frame angle about the parent's positive z-axis, radians |
| `c_linkage_dynamic_solve_joint`, `DYN_PRISMATIC_JOINT` | Signed attachment-origin separation along the parent joint frame's positive z-axis, model length units |
| `c_linkage_dynamic_solve` | Absolute moving-body world-z angle, radians; planar models only |

The joint-based entry point supports serial, planar, and spatial models, including
moving parents and arbitrarily oriented sliding axes. The library constructs the
driver internally; no constraint callback or axis argument is required.
Coincident joint frames define zero, not the initial configuration.

For example, after constructing a model whose second joint is prismatic and has
initial displacement `0.25`:

```c
static double displacement(double t, void *user_data)
{
    const double start = *(const double *)user_data;
    return start + 0.01 * (1.0 - cos(2.0*t)); /* zero initial rate */
}

/* In the calling function; settings and histories have been allocated. */
double start = 0.25;
int nbody = c_linkage_dynamic_body_count(model);
int nconstraint = c_linkage_dynamic_constraint_count(model) + 1;
c_linkage_dynamic_solve_joint(model, nbody, nconstraint, &settings, ntime, dt,
    gravity, body_force, body_torque, 2, displacement, &start,
    position, orientation, velocity, angular_velocity, multipliers,
    &converged, &iterations, &jacobian_singular, &completed_steps);
```

Include `<math.h>` for `cos`. Use one-based **joint** indices, not moving-body or
joint-variable indices: serial joint `i` is link `i`'s proximal joint; other models
use the original mechanism joint order. The joint's `actuated` flag is not required.
Only one input per solve is supported. Unsupported joint types, such as cylindrical
or spherical, and out-of-range indices are reported through `c_get_last_error`.

With a NULL callback, the selector is ignored and `nconstraint` must equal the
base count. With a non-NULL callback, allocate exactly one additional multiplier
row. The last row is torque for a revolute input or force for a prismatic input,
acting on the child relative to the parent toward increasing coordinate.
`c_linkage_dynamic_joint_reactions` excludes this actuator load.

The initial coordinate must match the callback at time zero within
`settings.tolerance`. The C linkage solver starts with the model's stored pose and
zero velocities, so choose a zero-initial-rate law for smooth startup. Initial
coordinate incompatibility is a deep Fortran input error and terminates the
process; C-level index/type/count errors return normally with an error code.
Angular residuals are wrapped modulo `2*pi`, permitting branch crossings and
multiple turns; resolve less than `pi` radians per time step.
The callback must be deterministic because it is evaluated repeatedly, including
for the final noncommitting look-ahead step. A dependent driver may overconstrain
a closed mechanism; inspect the convergence outputs rather than assuming success.

## Reference

The generated reference is organized by API area:

- [Public constants](./group__dynamics__constants.html)
- [Matrix and general kinematics](./group__dynamics__matrix.html)
- [Frequency response and system identification](./group__dynamics__frequency.html)
- [Quaternion operations](./group__dynamics__quaternion.html)
- [Geometry operations](./group__dynamics__geometry.html)
- [Serial linkage operations](./group__dynamics__serial.html)
- [Parallel and planar linkage operations](./group__dynamics__parallel.html)
- [Variational and linkage dynamics](./group__dynamics__variational.html)
- [Transfer functions and state-space models](./group__dynamics__state.html)
- [Structural analysis and line elements](./group__dynamics__structural.html)

The complete declaration reference is also available from the generated [Data Structures](annotated.html) and [Files](files.html) pages.

Routine-level documentation, including the purpose of every public routine, its
arguments, output arguments, and return value, is maintained in the [C API
routine reference](c-api-reference.md). The declarations in `dynamics.h` remain
the authoritative source for types, array dimensions, and calling conventions.

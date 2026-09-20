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
	dynamic-model handles. `c_linkage_dynamic_solve` then applies gravity,
	constant body loads, and optional prescribed planar motion.

Initialize `c_variational_integrator_settings` with
`c_default_variational_integrator_settings` before changing individual fields.
The direct and linkage solve routines write caller-owned histories using these
column-major shapes:

- position, velocity, and angular velocity: `3 × nbody × ntime`;
- orientation: `nbody × ntime`;
- multipliers: `nconstraint × ntime`.

Multiplier column `i` belongs to simulation point `i`. The final column is
computed with a noncommitting look-ahead step. For prescribed linkage motion,
the prescribed-angle multiplier is appended after the linkage constraints and
represents the required actuator torque.

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

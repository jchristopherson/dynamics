# DYNAMICS C API Routine Reference

The public C interface is declared in `src/c/dynamics.h`. Every routine is
listed in the generated Doxygen reference by API group. The declaration page
is the authoritative reference for the exact prototype; this page records the
meaning of the argument and return-value conventions used by every routine.

## Reading routine declarations

Each routine's declaration identifies every argument by name and C type. Input
arguments are marked `const` where the routine does not modify the referenced
storage. Non-const pointer and array arguments are caller-provided output or
in/out storage. Their required length or matrix shape is given by the related
size argument and the routine's API-group documentation.

For routines returning `void`, results are written to those non-const output
arguments. Routines returning `double`, `bool`, or `int` return the computed
scalar or status/count value shown in the declaration. Allocation routines
return an integer status and initialize the supplied object; a nonzero status
indicates allocation or construction failure. `c_alloc_dynamic_system_measurement_array`
returns the allocated measurement-array pointer, and linkage creation routines
return a `c_mechanism` handle. Release routines return `void` and accept the
object they own.

## Argument conventions

- `n`, `m`, `k`, `nfreq`, and similar integer arguments specify vector lengths,
  matrix dimensions, counts, or numbers of samples as named by the routine.
- `ld*` arguments are leading dimensions for column-major matrices. Matrix
  element `(i, j)` is stored at `j * ld + i`.
- Arrays marked `const` are read-only inputs. Arrays and structures without
  `const` are filled or updated by the routine and must have sufficient space.
- `c_iteration_behavior`, `c_regression_statistics`, and related structure
  pointers are output structures unless marked `const`.
- Mechanism link, joint, and frame indices are one-based. A `c_mechanism`
  returned by a creation routine must be released with `c_free_mechanism`.
- A `c_linkage_dynamic_model` returned by a creation routine must be released
  with `c_free_linkage_dynamic_model`.
- Callback arguments describe caller-supplied functions. Their callback typedef
  declarations in `dynamics.h` document the callback arguments and outputs.

## Variational history conventions

`c_variational_integrator_solve` and `c_linkage_dynamic_solve` use
column-major history arrays. Translational and angular histories have shape
`3 × nbody × ntime`, quaternion histories have shape `nbody × ntime`, and
multiplier histories have shape `nconstraint × ntime`.

The multiplier at column `i` is the discrete constraint multiplier associated
with the interval beginning at simulation point `i`. The final multiplier is
evaluated by a noncommitting look-ahead step so that forces and reactions are
defined at all returned simulation points.

For a prescribed planar-link motion, the total constraint count is the value
returned by `c_linkage_dynamic_constraint_count` plus one. The last multiplier
row is the required actuator torque. The prescribed-motion callback receives
simulation time and the caller's `user_data` pointer.

## Linkage reaction conventions

`c_linkage_dynamic_joint_reactions` accepts one state and its matching
multiplier column. It returns one `c_joint_reaction` per mechanism joint. Force
and moment vectors are expressed in world coordinates and act on the child
link. The equal-and-opposite wrench acts on the parent link. Multipliers for
planar-body restrictions and prescribed motion are not reported as joint
reactions.

## Complete routine list

The generated pages linked from the [C API guide](c-api.md) provide the complete
routine list and exact signatures, grouped as follows:

- Matrix and general kinematics
- Frequency response and system identification
- Quaternion operations
- Geometry operations
- Serial linkage operations
- Parallel and planar linkage operations
- Variational and linkage dynamics
- Transfer functions and state-space models

The routine names, associated argument names and types, and return types are
therefore documented together in the generated declaration reference rather
than duplicated here. This keeps the external C documentation synchronized
with `dynamics.h` when the interface evolves.

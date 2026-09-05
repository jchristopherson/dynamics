# DYNAMICS C API

The DYNAMICS C API is a C-compatible interface to the library. The public declarations are in [`dynamics.h`](../src/c/dynamics.h).

## Conventions

- Matrices use column-major storage, matching the Fortran API. For a matrix with `m` rows, `n` columns, and leading dimension `ld`, element `(i, j)` is stored at `j * ld + i` using zero-based C indices.
- Mechanism link, joint, and frame indices are one-based, matching the Fortran API.
- Memory returned by an allocation or creation routine is owned by the caller and must be released with its corresponding `c_free_*` routine.
- An optional pointer may be passed as `NULL` only where the function declaration or reference explicitly permits it.
- Solver callbacks use C calling conventions and receive array dimensions as explicit arguments.

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
[`examples/c_four_bar_example_1.c`](../examples/c_four_bar_example_1.c). It shows
how to allocate mechanism links, connect them with `c_joint` values, create a
planar linkage, solve its forward kinematics, and release all allocated data.

## Reference

The generated reference is organized by API area:

- [Public constants](./group__dynamics__constants.html)
- [Matrix and general kinematics](./group__dynamics__matrix.html)
- [Frequency response and system identification](./group__dynamics__frequency.html)
- [Quaternion operations](./group__dynamics__quaternion.html)
- [Geometry operations](./group__dynamics__geometry.html)
- [Serial linkage operations](./group__dynamics__serial.html)
- [Parallel and planar linkage operations](./group__dynamics__parallel.html)
- [Transfer functions and state-space models](./group__dynamics__state.html)

The complete declaration reference is also available from the generated [Data Structures](annotated.html) and [Files](files.html) pages.

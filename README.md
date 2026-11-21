# dynamics
A library of routines used for the analysis of dynamic systems.

## Status
[![CMake](https://github.com/jchristopherson/dynamics/actions/workflows/cmake.yml/badge.svg)](https://github.com/jchristopherson/dynamics/actions/workflows/cmake.yml)
[![Actions Status](https://github.com/jchristopherson/dynamics/workflows/fpm/badge.svg)](https://github.com/jchristopherson/dynamics/actions)

## Feature Matrix

| Category | Functionality | Status | Documentation | Examples |
|-----------|---------------|----------|---------------|----------|
| Frequency Response Analysis | SISO/MIMO FRFs, nonlinear frequency sweeps, modal response | ✅ | ✅ | ✅ |
| Modal Analysis | Modal parameter extraction, damping estimation, FRF fitting | ✅ | ✅ | ✅ |
| System Identification | Least-squares model fitting, regression statistics | ✅ | ✅ | ✅ |
| Controls & LTI Systems | State-space models, transfer functions, simulation tools | ✅ | ✅ | ✅ |
| Vibration Analysis | Q-factor, damping ratio, settling metrics, free-response characterization | ✅ | ✅ | ✅ |
| Kinematics | Denavit-Hartenberg models, forward and inverse kinematics | ✅ | ✅ | ✅ |
| Robotics Utilities | Jacobians, transformations, quaternion algebra | ✅ | ✅ | ✅ |
| Serial-Link Mechanisms | Revolute/prismatic linkages and motion analysis | ✅ | ✅ | ✅ |
| Parallel Mechanisms | Loop-closure constraints, mobility analysis, Jacobian partitioning | ✅ | ✅ | ✅ |
| Multibody Dynamics | Rigid-body variational integration | ✅ | ✅ | ✅ |
| Force Elements | Springs, dampers, torsional elements, custom force models | ✅ | ✅ | ✅ |
| Geometry | Points, lines, planes, Plücker lines, intersections, projections | ✅ | ✅ | ✅ |
| Structural Mechanics | Beam, truss, and shell finite elements | ✅ | ✅ | ✅ |
| Structural Dynamics | Generalized-α integration, dense and sparse assemblies | ✅ | ✅ | ✅ |
| Stability Analysis | Fixed-point stability classification | ✅ | ✅ | ✅ |
| C Interface | Native C API | ✅ | ✅ | ✅ |

## Documentation
- [Fortran API documentation](https://jchristopherson.github.io/dynamics/)
- [C API documentation](https://jchristopherson.github.io/dynamics/c-api/)

### Examples and Module Guides

See the subdirectory READMEs for detailed module descriptions, code examples, and plots where available.

- [C API](src/c/README.md) — C calling conventions and examples for kinematics and LTI simulation.
- [Controls](src/controls/README.md) — state-space models, transfer functions, and PID control.
- [Flexible Bodies](src/flexible_bodies/README.md) — flexible beams, material models, and mesh assembly.
- [Frequency Analysis](src/frequency/README.md) — modal analysis, frequency response, and nonlinear sweeps.
- [Rigid Bodies and Mechanisms](src/rigid_bodies/README.md) — serial and closed-loop kinematics and linkage dynamics.
- [Solvers](src/solvers/README.md) — linear structural time integration and constrained rigid-body integration.
- [Structural Analysis](src/structural/README.md) — beam, truss, shell, and discrete elements.
- [System Identification](src/system_id/README.md) — least-squares fitting of ODE models to measured data.
- [Vibrations and Stability](src/vibrations/README.md) — response metrics, fixed-point stability, and Poincaré sections.

## Building

### Prerequisites
- A Fortran compiler with Fortran 2018 support.
- CMake 3.24+ for the CMake workflow.
- fpm for the FPM workflow.
- BLAS and LAPACK available to the linker (required by the FPM build configuration).

### Build With CMake
Configure and build from the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

Useful CMake options:
- `-DBUILD_TESTING=ON` to build test targets.
- `-DBUILD_DYNAMICS_EXAMPLES=ON` to build example programs.
- `-DBUILD_DYNAMICS_C_INTERFACE=ON` to include the C interface.
- `-DBUILD_SHARED_LIBS=ON` to build shared libraries.

Run tests (if enabled):

```sh
ctest --test-dir build --output-on-failure
```

Install the library:

```sh
cmake --install build --prefix <install-prefix>
```

### Build With FPM
From the repository root:

```sh
fpm build --profile release
```

Run tests:

```sh
fpm test
```

Install:

```sh
fpm install --prefix <install-prefix>
```

FPM resolves the package dependencies declared in `fpm.toml` automatically.

## Contributing
Contributors should follow the [style guide](STYLE_GUIDE.md) for Fortran source, tests, and documentation comments.

## Quick Start

### Use From FPM
Add `dynamics` to your `fpm.toml` dependencies:

```toml
[dependencies]
dynamics = { git = "https://github.com/jchristopherson/dynamics.git" }
```

Then build and run your project:

```sh
fpm build
fpm run
```

### Use From CMake
If `dynamics` is installed and discoverable via `CMAKE_PREFIX_PATH`, link it as a package:

```cmake
find_package(dynamics REQUIRED)
target_link_libraries(your_target PRIVATE dynamics::dynamics)
```

If you prefer vendoring it directly, add it as a subdirectory and link the target:

```cmake
add_subdirectory(path/to/dynamics)
target_link_libraries(your_target PRIVATE dynamics)
```


## References
1. J. D. Hartog, "Mechanical Vibrations," New York: Dover Publications, Inc., 1985.
2. S. S. Rau, "Mechanical Vibrations," 3rd ed., Reading, MA: Addison-Wesley Publishing Co., 1995.
3. R. N. Jazar, "Advanced Vibrations," 2nd ed., New York: Springer, 2022.
4. W. T. Thomson, "Theory of Vibration with Applications," 4th ed., New York: Springer, 1993.
5. A. H. Nayfeh and B. Balachandran, "Applied Nonlinear Dynamics. Analytical, Computational, and Experimental Methods," New York: John WIley & Sons, Inc., 1995.
6. L. Meirovitch, "Fundamentals of Vibrations," Long Grove, IL: Waveland Press, Inc., 2001.
7. R. N. Jazar, "Theory of Applied Robotics, Kinematics, Dynamics, and Control," New York: Springer, 2007.
8. A. H. Nayfeh, "Introduction to Perturbation Techniques," New York: John Wiley & Sons, Inc., 1993.
9. Jolicoeur, M.P., Roumy, J.G., Vanreusel, S., Dionne, D., Douville, H., Boulet, B., Michalska, H., Masson, P., & Berry, A. (2005). "Reduction of structure-borne noise in automobiles by multivariable feedback." 1397 - 1402. 10.1109/CCA.2005.1507327. 
10. Brunton, Steven & Proctor, Joshua & Kutz, J.. (2015). "Discovering governing equations from data: Sparse identification of nonlinear dynamical systems." Proceedings of the National Academy of Sciences. 113. 3932–3937. 10.1073/pnas.1517384113. 
11. Brüdigam, Jan & Sosnowski, Stefan & Manchester, Zac & Hirche, Sandra. (2023). Variational integrators and graph-based solvers for multibody dynamics in maximal coordinates. Multibody System Dynamics. 61. 1-34. 10.1007/s11044-023-09949-x. 
12. J.-L. Batoz, K.-J. Bathe, and L.-W. Ho, "A study of three-node triangular plate bending elements," International Journal for Numerical Methods in Engineering, vol. 15, no. 12, pp. 1771-1812, 1980. 10.1002/nme.1620151205.
13. E. N. Dvorkin and K.-J. Bathe, "A continuum mechanics based four-node shell element for general non-linear analysis," Engineering Computations, vol. 1, no. 1, pp. 77-88, 1984. 10.1108/eb023562.
14. T. J. R. Hughes and F. Brezzi, "On drilling degrees of freedom," Computer Methods in Applied Mechanics and Engineering, vol. 72, no. 1, pp. 105-121, 1989. 10.1016/0045-7825(89)90124-2.
15. R. D. Cook, D. S. Malkus, M. E. Plesha, and R. J. Witt, "Concepts and Applications of Finite Element Analysis," 4th ed., New York: John Wiley & Sons, Inc., 2002.

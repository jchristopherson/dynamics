## Structural Analysis Modules

The structural finite-element functionality is organized into four element modules and a shared foundation. The modules can be used directly, or through the `dynamics` umbrella module used by the examples.

* **`dynamics_structural`** provides the common finite-element types and operations: `node`, isotropic `material`, abstract `element` and `line_element` base types, element shape-function utilities, and global static/dynamic matrix assembly. It also provides boundary-condition helpers, restoration of constrained degrees of freedom, static-system solving, and nodally averaged strain and stress. Assembly and several utilities support dense and CSR sparse matrices. The integration-rule constants select one- through four-point integration.
* **`dynamics_beam_elements`** provides two-node Bernoulli-Euler beam elements in 2D and 3D. They include bending response as well as stiffness, mass, strain, stress, shear-force, and bending-moment calculations. The 3D element uses an orientation point to define its local axes.
* **`dynamics_truss_elements`** provides 2D and 3D two-node pin-jointed bars. These elements carry axial load only, with translational degrees of freedom at each node; they do not model bending or nodal rotations.
* **`dynamics_discrete_elements`** provides lumped-parameter axial springs and viscous dampers between two nodes, and translational point masses attached to one node, in 2D and 3D. Its assembly routines combine these contributions into mass, damping, and stiffness matrices. These elements use translational degrees of freedom, mapping to the leading translational components when attached to nodes with additional degrees of freedom.
* **`dynamics_shell_elements`** provides flat 3D triangular and quadrilateral shell elements with six degrees of freedom per node (three translations and three rotations). The triangle is intended for thin shells and neglects transverse shear deformation; the quadrilateral uses an MITC4 formulation and supports thin and thick shells. Node ordering defines the element normal.

These element modules depend on the shared `node`, `material`, and assembly interfaces from `dynamics_structural`. Typically, create nodes and elements, assemble the global matrices, remove constrained degrees of freedom with `apply_boundary_conditions`, then pass the reduced system to a static, modal, frequency-response, or time-integration analysis.

## Structural Example

The [`structural_example`](examples/structural_example.f90) analyzes a two-dimensional aluminum frame shaped as a triangular truss with a center member. Each of its five members is divided into 50 Euler-Bernoulli beam elements, giving 250 elements and 249 nodes. Each node has two translational degrees of freedom and one in-plane rotation. The lower-left node is fixed in all three degrees of freedom, the lower-right node is supported vertically, and a 10 kN downward load is applied at the apex.

The example assembles global stiffness and mass matrices, solves the constrained static system, recovers support reactions, and computes the six lowest modes. It then forms a Rayleigh-damped frequency response from 10 Hz to 1 kHz, exciting the apex vertically. The central analysis calls are:

```fortran
call assemble_dynamic_system(gdof, beams, nodes, M, K)

allocate(F(size(K, 1)), source = 0.0d0)
F(5) = -applied_force
bc = [1, 2, 3, 8]
Fbc = apply_boundary_conditions(bc, F)
Kbc = apply_boundary_conditions(bc, K)
Mbc = apply_boundary_conditions(bc, M)

ubc = solve_static_system(Kbc, Fbc)
u = restore_constrained_values(bc, ubc)
call modal_response(Mbc, Kbc, n_plot_modes, freqs, shapesbc)

excitefcn => modal_frf_forcing_term
frsp = frequency_response(Mbc, Kbc, alpha, beta, n_plot_modes, nfreq, &
    minfreq, maxfreq, excitefcn)
```

The example uses SI units: the geometry is in metres, material properties are in pascals and kg/m^3, and the applied force is in newtons. Each node's global degrees of freedom are ordered as \(x\)-translation, \(y\)-translation, and in-plane rotation. Thus, the downward force at node 2 is assigned to global degree of freedom 5, and `bc = [1, 2, 3, 8]` fixes node 1 and the vertical degree of freedom at node 3. `apply_boundary_conditions` removes these indices from the system; `restore_constrained_values` puts zero values back into the reduced static solution and mode shapes. The reactions are then recovered from the original, unconstrained stiffness matrix.

The static plot overlays the undeformed frame with the deformed shape, magnified for visibility. The text output below reports the first four nodes' displacement components multiplied by 1000 and the summed support reactions.

![Static deformation of the five-member frame; deformation is magnified for visibility](../../images/structural_example_static_deformation.png?raw=true)

```txt
NODAL DISPLACEMENTS (x 1000)
Node       X              Y                Theta
1       0.000E+0        0.000E+0        0.000E+0
2       -41.247E-3      -229.425E-3     -60.831E-3
3       -81.263E-3      0.000E+0        692.776E-3
4       -39.967E-3      -228.573E-3     -65.412E-3

REACTION LOADS
Node      FX              FY               MZ
1       9.234E-9        5.209E+3        208.794E+0
3       600.829E-12     4.791E+3        1.364E-12
SUM     9.834E-9        10.000E+3       208.794E+0

MODAL RESPONSE
Mode 1: 379.810 Hz
Mode 2: 553.835 Hz
Mode 3: 670.329 Hz
Mode 4: 804.472 Hz
Mode 5: 949.697 Hz
Mode 6: 1283.488 Hz
```

The first six natural frequencies and corresponding mode shapes are:

| Mode 1: 379.810 Hz | Mode 2: 553.835 Hz | Mode 3: 670.329 Hz |
|:--:|:--:|:--:|
| ![Mode 1](../../images/structural_example_mode_1.png?raw=true) | ![Mode 2](../../images/structural_example_mode_2.png?raw=true) | ![Mode 3](../../images/structural_example_mode_3.png?raw=true) |
| **Mode 4: 804.472 Hz** | **Mode 5: 949.697 Hz** | **Mode 6: 1283.488 Hz** |
| ![Mode 4](../../images/structural_example_mode_4.png?raw=true) | ![Mode 5](../../images/structural_example_mode_5.png?raw=true) | ![Mode 6](../../images/structural_example_mode_6.png?raw=true) |

For the harmonic analysis, the example uses proportional damping, $C=\alpha M+\beta K$, with $\alpha=10^{-3}$ and $\beta=2\times10^{-6}$. The following plot shows the normalized magnitude and phase of the apex vertical response (degree of freedom 5) under a 10 kN harmonic force.

![Frequency response of the apex vertical degree of freedom from 10 Hz to 1 kHz](../../images/structural_example_frf.png?raw=true)

The modal solver returns angular frequencies, which the example converts from rad/s to Hz for reporting and plotting. The frequency-response calculation sweeps 1000 points from 10 Hz to 1 kHz and uses the first six modes. The deformed-shape plot applies an additional visualization scale; it is not drawn at physical displacement scale.

Build with `BUILD_DYNAMICS_EXAMPLES=ON`, then build and run the `structural_example` target.


## Harmonic Truss Example

The [`harmonic_truss_example`](examples/harmonic_truss_example.f90) models a four-node, five-bar pin-jointed truss with a pinned support and a roller. It assembles axial stiffness and consistent translational mass matrices, applies the support constraints, and advances the response to a vertical 20 Hz load at the apex with `dense_generalized_alpha_integrator`. Fplot saves vertical apex and midspan motion and horizontal midspan and roller motion to `harmonic_truss_response.png` in the process's working directory. Build the `harmonic_truss_example` CMake target with examples enabled, then run it from the build's examples directory.

After assembling and reducing the mass and stiffness matrices, the time-stepping core is:

```fortran
call integrator%initialize(reduced_mass, damping, reduced_stiffness, rho_infinity = 0.7d0)
do time_index = 1, nsteps + 1
    time(time_index) = (time_index - 1) * dt
    apex_motion(time_index) = 1.0d3 * displacement(apex_y)
    middle_motion(time_index) = 1.0d3 * displacement(middle_y)
    if (time_index > nsteps) exit

    force_next = 0.0d0
    force_next(apex_y) = -force_amplitude * &
        sin(2.0d0 * pi * frequency * real(time_index, real64) * dt)
    call integrator%step(force_current, force_next, dt, &
        displacement, velocity, acceleration)
    force_current = force_next
end do
```

This example uses 2 translational degrees of freedom per node and axial-only bars, so the five bars do not carry bending moments. The restrained global degrees of freedom are `[1, 2, 6]` (both translations at the pinned node and the vertical translation at the roller). Since constrained entries are removed, the example adjusts the global indices when locating the apex and sensor degrees of freedom in the reduced vectors. The integrator advances 2000 steps of 0.2 ms, and the plotted displacements are converted from metres to millimetres.

![Forced response of the simple truss example](../../images/harmonic_truss_example.png?raw=true)

## Shell Element Example

The [`shell_element_example`](examples/shell_element_example.f90) computes the first six modes of a flat 0.75 m by 1.25 m aluminum plate, 10 mm thick, clamped along its two long edges ($x=0$ and $x=0.75$ m). The plate is meshed with a 50-by-50 grid of nodes and 2401 `rectangular_shell_element` objects. Each node has six degrees of freedom: three translations and three rotations. Fplot draws each mode as a surface of the out-of-plane ($z$) displacement over the mesh.

The core of the example builds the mesh from a `meshgrid`, assembles sparse mass and stiffness matrices, clamps the two edges, and solves for the modes:

```fortran
! Create the mesh
xc = linspace(0.0d0, width, ndiv)
yc = linspace(0.0d0, length, ndiv)
xy = meshgrid(xc, yc)
kk = 0
do j = 1, ndiv
    do i = 1, ndiv
        kk = kk + 1
        nodes(kk) = node(kk, dof_per_node, xy(i,j,1), xy(i,j,2), 0.0d0)
    end do
end do

! Nodes are listed counter-clockwise to define the element normal
elements(kk) = rectangular_shell_element(mat, thickness, &
    nodes(n1), nodes(n2), nodes(n3), nodes(n4))

! Assemble, constrain, and solve
call assemble_dynamic_system(global_dof_count, elements, nodes, M, K)
Mbc = apply_boundary_conditions(bc, M)
Kbc = apply_boundary_conditions(bc, K)
call modal_response(Mbc, Kbc, mode_count, freqs, shapesbc)
freqs = freqs / (2.0d0 * pi)

! Restore the constrained DOF and plot the z-displacement of each mode
do i = 1, mode_count
    shapes(:,i) = restore_constrained_values(bc, shapesbc(:,i))
    call extract_nodal_displacements(shapes(:,i), u)
    call plot_mode_shape(i, freqs(i), xy(:,:,1), xy(:,:,2), u)
end do
```

The element nodes are supplied counter-clockwise so their ordering defines a consistent local normal. The two clamped edges constrain all six degrees of freedom at every node on those edges. Modal analysis is performed on the reduced sparse matrices; constrained zeros are restored to each mode shape before the nodal displacement extractor selects the out-of-plane component for plotting. The reported frequencies are converted from rad/s to Hz.

```txt
MODAL RESPONSE:
Mode 1: 97.801 Hz
Mode 2: 104.521 Hz
Mode 3: 130.581 Hz
Mode 4: 183.377 Hz
Mode 5: 268.814 Hz
Mode 6: 269.948 Hz
```

The first six natural frequencies and corresponding mode shapes are:

| Mode 1: 97.801 Hz | Mode 2: 104.521 Hz | Mode 3: 130.581 Hz |
|:--:|:--:|:--:|
| ![Mode 1](../../images/shell_example_mode_1.png?raw=true) | ![Mode 2](../../images/shell_example_mode_2.png?raw=true) | ![Mode 3](../../images/shell_example_mode_3.png?raw=true) |
| **Mode 4: 183.377 Hz** | **Mode 5: 268.814 Hz** | **Mode 6: 269.948 Hz** |
| ![Mode 4](../../images/shell_example_mode_4.png?raw=true) | ![Mode 5](../../images/shell_example_mode_5.png?raw=true) | ![Mode 6](../../images/shell_example_mode_6.png?raw=true) |
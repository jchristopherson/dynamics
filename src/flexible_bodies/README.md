## Flexible-Body Modules

* **`dynamics_flexible_body_types`** defines `fe_node` reference coordinates and the `degree_of_freedom` bookkeeping type.
* **`dynamics_flexible_element`** defines the common abstract element contract for degree-of-freedom mapping, mass, internal force, and tangent stiffness.
* **`dynamics_finite_element`** extends that contract with mesh connectivity, node count, degrees of freedom per node, and spatial dimensionality.
* **`dynamics_flexible_materials`** defines strain-energy constitutive interfaces and the built-in linear-elastic isotropic material. Material routines provide derivatives through real, dual, and hyper-dual evaluations for force and tangent calculations.
* **`dynamics_flexible_beam_elements`** provides energy-based spatial beam elements, including the two-node Euler-Bernoulli beam described below.
* **`dynamics_flexible_assembly`** provides heterogeneous finite-element meshes, connectivity-based equation numbering, and dense or CSR assembly of mass, tangent stiffness, and internal force.
* **`dynamics_boundary_conditions`** applies force and displacement boundary conditions to a mesh's assembled equations, dense or CSR, and restores prescribed displacements after the reduced solve.

## Flexible beams

`euler_bernoulli_beam_element` is a two-node spatial, shear-free beam with no
prescribed section shape. Store reference node coordinates along the centroidal
axis in a caller-owned `fe_node` array. Set each beam's `node_indices` to the two
one-based positions in that array, not the nodes' `id` values. Pass the same
array to `mass_matrix`, `strain_energy`, `internal_force_vector`, and
`stiffness_matrix`. Beams can share an index and will read any coordinate
updates on the next evaluation. Remap connectivity if the array is reordered.
Set the beam's material pointer and reference local-y direction along a principal
axis. Supply finite, positive `area`, `j` (Saint-Venant torsion constant),
`iyy`, and `izz` (second moments of area about the local principal axes).
Assign `equation_indices` only if you need them for manual assembly. Each node has
global `[ux, uy, uz, rx, ry, rz]` DOFs; rotations are rotation vectors in
radians. The material target must remain alive while the beam uses it.

The beam integrates section energy at four Gauss points along the length.
It evaluates the material's `plane_stress_energy` at the centroidal
Green-Lagrange strain and multiplies it by `area`. The material supplies
`bending_energy(ky, kz, iyy, izz)` and `torsional_energy(twist, j)` per unit
reference length. These laws may be nonlinear, but axial, bending, and
torsional energies are separate; there is no through-section fiber integration.
The consistent reference mass uses `iyy + izz` for polar rotary inertia, not `j`.
Internal force is the positive energy gradient and stiffness its Hessian;
the beam's AD workspace is reused and must not be accessed concurrently.

### Flexible mesh assembly

See the [Flexible Truss Example](#flexible-truss-example) for a complete
sparse static and modal analysis.

Populate a `finite_element_mesh` with its `nodes` array and allocate each
`elements(i)%element` from a concrete finite element. Elements reference nodes
by array position. `number_dofs(offsets)` returns the global equation ranges:
node `i` uses `offsets(i)+1` through `offsets(i+1)`, and unused nodes take no
equations. Concrete flexible finite elements must implement `dof_components()`
with distinct positive labels in their local node-block order. Use `DOF_UX`
through `DOF_RZ` for the standard mechanical components. If you implement your
own elements, add this binding when updating to the new interface. The mesh
rejects elements that disagree on DOF counts, labels, or label order at a shared
node. Assembly does not modify element equation indices.

`mesh%resolve_dofs(node_positions, components, dofs[, ndof])` returns
`degree_of_freedom` records in request order. Inputs use node-table positions;
each result reports the node's `fe_node%id`, its physical component label, and
its current global equation number. Resolve again if the mesh layout changes.

Call `assemble_mass(matrix)` or `assemble_stiffness(q, matrix)` with **either**
an allocatable `real(real64)` rank-two matrix **or** a `csr_matrix` output.
The CSR path assembles directly into sparse storage. Use `assemble_internal_force(q, f)`
for the global internal-force vector. Both displacement `q` and force `f` use
the mesh equation ordering. These outputs are assembled before applying
boundary conditions.

With OpenMP, pass `nthreads=4` (for example) to `assemble_stiffness` or
`assemble_internal_force` to evaluate elements concurrently. The default is
one thread, and builds without OpenMP run serially. Mass assembly is always
serial. Contributions are added in element order, even when evaluation is
parallel. Only opt in if your element and material routines are safe to call
concurrently on distinct elements. The built-in Euler–Bernoulli beam has a
workspace per element and can share a read-only built-in material. Do not
assemble the same mesh concurrently from multiple callers.

For example, given a configured `beam` and its `node_table`, assemble the
reference mass and the tangent at zero displacement:

```fortran
! Declare alongside beam and node_table.
type(finite_element_mesh) :: mesh
integer, allocatable, dimension(:) :: offsets
real(real64), allocatable, dimension(:) :: q
real(real64), allocatable, dimension(:) :: force
type(csr_matrix) :: mass
type(csr_matrix) :: stiffness

mesh%nodes = node_table
allocate(mesh%elements(1))
allocate(mesh%elements(1)%element, source=beam)
call mesh%number_dofs(offsets)
allocate(q(offsets(size(offsets))), source=0.0_real64)
call mesh%assemble_mass(mass)
call mesh%assemble_stiffness(q, stiffness)
call mesh%assemble_internal_force(q, force)
```

### Mesh boundary conditions

`dynamics_boundary_conditions` addresses a DOF with a `nodal_condition`, which
holds a one-based node position in `mesh%nodes`, a physical component label
(e.g., `DOF_UY`, not the component's local slot), and a `value`.
`resolve_boundary_conditions(mesh, conditions, dofs[, ndof])` returns the
global equation indices without applying the conditions. For homogeneous
supports, use those indices with `apply_boundary_conditions` to reduce multiple
matrices, as in the flexible truss example. Force conditions accumulate into
the global force vector, which is allocated as zeros if unallocated:
`apply_force_boundary_conditions(mesh, conditions, f)`.

Displacement conditions are eliminated symmetrically from the assembled system
with `apply_displacement_boundary_conditions(mesh, conditions, k, f)`, where `k`
is dense or `csr_matrix`. The matrix and force are reduced to the free DOFs, and
the force becomes `f - K*u_p`. Duplicate displacement DOFs are rejected. After
solving the reduced system, `restore_displacement_boundary_conditions(mesh,
conditions, uf)` returns the full displacement vector with the prescribed values.
Invalid node or component indices, and size mismatches, cause `error stop`.

For dense output, declare `mass` and `stiffness` as allocatable rank-two
`real(real64)` arrays instead; the calls are unchanged. Allocating from `beam`
copies the element but not its material target, which must remain alive while
the mesh uses it.

## Flexible Truss Example

The [`flexible_truss_example`](../../examples/flexible_truss_example.f90)
reproduces the five-member geometry of the
[structural example](../structural/README.md#structural-example) using spatial
Euler–Bernoulli beams. Each member is divided into 50 elements, giving 250
elements and 249 nodes. All six DOFs are fixed at the lower-left support
(node 1); only vertical translation (`DOF_UY`) is restrained at the
lower-right support (node 3). The mesh has 1,487 free DOFs.

![Undeformed five-member flexible truss mesh with markers at the beam nodes](../../images/flexible_truss_example_mesh.png?raw=true)

The example assembles the reference mass and undeformed tangent stiffness in
CSR format. The linearized static solution for a 10 kN downward load at the
apex gives an apex displacement of **0.229 mm downward**. The XY plot overlays
the dotted reference outline with the displacement-colored shape, magnified
100 times for visibility.

![Static deformation of the flexible truss under a 10 kN apex load, magnified 100 times](../../images/flexible_truss_example_static_deformation.png?raw=true)

The example also solves nonlinear equilibrium under the same load, using the
linear displacement as its initial guess. It passes the **full** load and
displacement vectors and the original nodal support conditions to the solver;
the solver applies the constraints itself. It reports the nonlinear apex
displacement, Newton iterations, and free residual norm, and displays a separate
deformation plot. The built-in material is linear elastic; the nonlinearity
here comes from the beam's finite-strain kinematics. For this load, the
nonlinear apex displacement is approximately **0.229436 mm downward**, close
to the linear result; the example converged in two Newton iterations. It uses
CSR tangents and ILUT-preconditioned GMRES for the nonlinear corrections.

```fortran
nonlinear_u = u  ! Full displacement recovered from the linear solution.
static_settings%sparse = .true.
call solve_flexible_static(mesh, force, supports, nonlinear_u, static_settings, static_info)
if (.not.static_info%converged) error stop "Nonlinear static equilibrium failed"
```

The six lowest vibration modes use the same undeformed tangent. The 3D
plots rotate the display coordinates by +90 degrees about X so global Y
points up and the depth axis is labeled `-Z`. Each mode is scaled
for display; the plotted amplitudes are not physical displacements.

| Mode 1: 33.173 Hz | Mode 2: 153.396 Hz | Mode 3: 251.955 Hz |
|:--:|:--:|:--:|
| ![Mode 1](../../images/flexible_truss_example_mode_1.png?raw=true) | ![Mode 2](../../images/flexible_truss_example_mode_2.png?raw=true) | ![Mode 3](../../images/flexible_truss_example_mode_3.png?raw=true) |
| **Mode 4: 379.459 Hz** | **Mode 5: 552.349 Hz** | **Mode 6: 554.742 Hz** |
| ![Mode 4](../../images/flexible_truss_example_mode_4.png?raw=true) | ![Mode 5](../../images/flexible_truss_example_mode_5.png?raw=true) | ![Mode 6](../../images/flexible_truss_example_mode_6.png?raw=true) |

Build the `flexible_truss_example` target with `BUILD_DYNAMICS_EXAMPLES=ON`
and run it with Gnuplot on `PATH`.
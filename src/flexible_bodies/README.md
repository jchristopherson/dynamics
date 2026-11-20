## Flexible-Body Modules

* **`dynamics_flexible_body_types`** defines `fe_node` reference coordinates and the `degree_of_freedom` bookkeeping type.
* **`dynamics_flexible_element`** defines the common abstract element contract for degree-of-freedom mapping, mass, internal force, and tangent stiffness.
* **`dynamics_finite_element`** extends that contract with mesh connectivity, node count, degrees of freedom per node, and spatial dimensionality.
* **`dynamics_flexible_materials`** defines strain-energy constitutive interfaces and the built-in linear-elastic isotropic material. Material routines provide derivatives through real, dual, and hyper-dual evaluations for force and tangent calculations.
* **`dynamics_flexible_beam_elements`** provides energy-based spatial beam elements, including the two-node Euler-Bernoulli beam described below.
* **`dynamics_flexible_assembly`** provides heterogeneous finite-element meshes, connectivity-based equation numbering, and dense or CSR assembly of mass, tangent stiffness, and internal force.

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

Populate a `finite_element_mesh` with its `nodes` array and allocate each
`elements(i)%element` from a concrete finite element. Elements reference nodes
by array position. `number_dofs(offsets)` returns the global equation ranges:
node `i` uses `offsets(i)+1` through `offsets(i+1)`, and unused nodes take no
equations. Shared nodes must agree on DOF count and component order. Mesh
assembly derives its own numbering without modifying element equation indices.

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

For dense output, declare `mass` and `stiffness` as allocatable rank-two
`real(real64)` arrays instead; the calls are unchanged. Allocating from `beam`
copies the element but not its material target, which must remain alive while
the mesh uses it.
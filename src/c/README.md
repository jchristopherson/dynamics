## C API

The public C interface is declared in [`dynamics.h`](dynamics.h). It exposes core DYNAMICS functionality to C programs, including:

* Matrix operations, coordinate transforms, Denavit-Hartenberg tables, and forward/inverse kinematics.
* Serial, planar, and spatial linkage modeling and analysis.
* Quaternion, geometry, and structural-analysis operations.
* Frequency response, system identification, stability, and vibration analysis.
* State-space and transfer-function models, plus time integration for linear systems.
* [Discrete-time controls](../../c_docs/c-api.md#discrete-time-controls), including sampled-model conversion, one-step evaluation, and integral state-feedback tuning.
* Variational integration of constrained rigid-body systems and linkage dynamics.

The examples below demonstrate three different parts of the API: constructing a Denavit-Hartenberg table, solving the kinematics of a closed-loop linkage, and integrating a continuous-time state-space model.

## Conventions and ownership

* Matrices use column-major storage, matching Fortran. For a matrix with `m` rows and leading dimension `ld`, element `(i, j)` uses zero-based C indices and is stored at `j * ld + i`.
* Mechanism link, joint, and frame indices are one-based, even though C arrays are zero-based.
* Allocation and creation routines return caller-owned storage or handles. Release each object with its matching `c_free_*` routine.
* Callback routines receive the `user_data` pointer passed by the caller, so application state can be carried without globals.
* C-level input errors can be checked with `c_get_last_error` and `c_get_last_error_message`. The C header notes that errors raised deeper in the Fortran library may terminate the process.

For detailed prototypes and array dimensions, see the [C API guide](../../c_docs/c-api.md), the [routine reference](../../c_docs/c-api-reference.md), and the declarations in [`dynamics.h`](dynamics.h).

## Denavit-Hartenberg table example

[`c_dh_table_example_1`](../../examples/c_dh_table_example_1.c) defines coordinate systems for a three-revolute-joint manipulator, builds a DH table from those frames, and evaluates the resulting end-to-end transformation. Once `csys` has been populated, the central operations are:

```c
/* Inside main, with csys[0..3] already defined. */
c_dh_table tbl;
double T[16];

/* Four coordinate systems describe the base and three link frames. */
c_build_dh_table(4, csys, &tbl);

/* T is a 4-by-4 homogeneous transform with leading dimension 4. */
c_dh_forward_kinematics_table(&tbl, T, 4);
printf("\nLinkage Transformation Matrix:\n");
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\n",
    T[0], T[4], T[8], T[12]);
printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\n",
    T[1], T[5], T[9], T[13]);

c_free_dh_table(&tbl);
```

In the complete example, each frame is assembled from the joint rotation and translated link offset. `c_build_dh_table` extracts the three parameter sets; `c_dh_forward_kinematics_table` composes their transforms. The transform matrix is column-major, so its translation occupies entries `T[12]`, `T[13]`, and `T[14]`.

## Four-bar linkage example

[`c_four_bar_example_1`](../../examples/c_four_bar_example_1.c) builds a planar four-bar linkage with four revolute joints. The example creates a starting configuration to select the desired assembly branch, sweeps the crank angle, and queries the coupler pose and Jacobian.

After defining each link's two joint frames, the linkage and a crank position can be set up and evaluated as follows:

```c
/* Inside main, with links[4], joints[4], and tool initialized. */
c_mechanism linkage;
c_iteration_behavior behavior;
double config[4] = {0.0, 0.5 * PI, -0.5 * PI, 0.0};
double theta = 0.5 * PI;
double T[16];
double jacobian[3];

/* link and joint indices in c_joint are one-based. */
linkage = c_create_planar_linkage(4, links, 4, joints, 1, 3, tool, 4);
c_mechanism_set_configuration(linkage, 4, config);
c_mechanism_forward_kinematics(linkage, 1, &theta, T, 4, &behavior);
c_mechanism_jacobian(linkage, 1, &theta, jacobian, 3);

c_free_mechanism(linkage);
```

The starting configuration matters because a closed-loop mechanism can have multiple valid assembly modes. The full example also reports the mechanism's loop count, constraints, and degrees of freedom, checks the loop-closure residual, and frees every allocated link with `c_free_mechanism_link`.

## LTI example

[`c_lti_example_1`](../../examples/c_lti_example_1.c) models a single-degree-of-freedom mass-spring-damper system, applies a sinusoidal force through a callback, and writes the computed response at the requested time samples.

```c
#include <math.h>
#include <dynamics.h>

static void excitation(int n, double t, double *u, void *user_data)
{
    u[0] = 1.0e3 * sin(10.0 * t);
}

int main(void)
{
    const double m = 1.0, b = 5.0, k = 2.0e3, tmax = 1.0;
    const int n = 100;
    double ic[2] = {0.0, 0.0}, t[100], y[100];
    double dt = tmax / (n - 1.0);
    c_state_space_model mdl;
    int i;

    c_create_state_space_model(1, 1, &m, 1, &b, 1, &k, 1, &mdl);
    for (i = 0; i < n; ++i) t[i] = i * dt;
    c_lti_solve(&mdl, excitation, n, t, 2, ic,
        DYN_RUNGE_KUTTA_45, 1, y, n, NULL);
    c_free_state_space_model(&mdl);
    return 0;
}
```

`c_create_state_space_model` constructs the first-order state-space realization of the second-order mechanical model. `c_lti_solve` advances that model using the selected integrator; the callback's `user_data` is `NULL` here because the input does not require extra parameters. See the full source for the time-grid construction and output loop.

## Building

Enable both the C interface and examples when configuring the project:

```sh
cmake -S . -B build -DBUILD_DYNAMICS_EXAMPLES=ON -DBUILD_DYNAMICS_C_INTERFACE=ON
cmake --build build --target c_dh_table_example_1 c_four_bar_example_1 c_lti_example_1
```

For an installed CMake consumer, include `dynamics.h` and link against the exported `dynamics::dynamics` target.
#include "dynamics_structural_tests.h"
#include "dynamics_c_test_helper.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static c_node make_node(int index, int dof, double x, double y, double z)
{
    c_node n;
    n.index = index;
    n.dof = dof;
    n.x = x;
    n.y = y;
    n.z = z;
    return n;
}

static void outer_product(int n, double scale, const double *b, double *x,
    int ldx)
{
    int i, j;
    for (j = 0; j < n; ++j)
        for (i = 0; i < n; ++i)
            x[INDEX(i, j, ldx)] = scale * b[i] * b[j];
}

bool c_test_truss_elements()
{
    // Local Variables
    const double tol = 1.0e-12;
    bool rst;
    double b2[4] = {-0.6, -0.8, 0.6, 0.8};
    double b3[6] = {-3.0 / 13.0, -4.0 / 13.0, -12.0 / 13.0, 3.0 / 13.0,
        4.0 / 13.0, 12.0 / 13.0};
    double u2[4] = {0.0, 0.0, 0.006, 0.008};
    double k2[16], m2[16], kans2[16], k3[36], kans3[36];
    c_material mat;
    c_node nodes2[2], nodes3[2];
    c_truss_element_2d bar2;
    c_truss_element_3d bar3;

    // Initialization
    rst = true;
    mat.modulus = 100.0;
    mat.poissons_ratio = 0.3;
    mat.density = 6.0;

    // 2D truss
    nodes2[0] = make_node(1, 2, 0.0, 0.0, 0.0);
    nodes2[1] = make_node(2, 2, 3.0, 4.0, 0.0);
    bar2.material = mat;
    bar2.area = 0.02;
    bar2.node_1 = nodes2[0];
    bar2.node_2 = nodes2[1];
    outer_product(4, 0.4, b2, kans2, 4);
    c_assemble_dynamic_system_truss_2d(4, 1, &bar2, 2, nodes2, 0, m2, 4, k2, 4);
    if (!compare_matrices(4, 4, k2, 4, kans2, 4, tol) ||
        fabs(m2[INDEX(0, 0, 4)] - 0.2) > tol ||
        fabs(m2[INDEX(0, 2, 4)] - 0.1) > tol ||
        fabs(c_truss_element_2d_length(&bar2) - 5.0) > tol ||
        fabs(c_truss_element_2d_strain(&bar2, u2) - 0.002) > tol ||
        fabs(c_truss_element_2d_axial_force(&bar2, u2) - 0.004) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_truss_elements - 2D\n");
    }

    // 3D truss
    nodes3[0] = make_node(1, 3, 0.0, 0.0, 0.0);
    nodes3[1] = make_node(2, 3, 3.0, 4.0, 12.0);
    bar3.material = mat;
    bar3.area = 0.02;
    bar3.node_1 = nodes3[0];
    bar3.node_2 = nodes3[1];
    outer_product(6, 2.0 / 13.0, b3, kans3, 6);
    c_assemble_static_system_truss_3d(6, 1, &bar3, 2, nodes3, k3, 6);
    if (!compare_matrices(6, 6, k3, 6, kans3, 6, tol) ||
        fabs(c_truss_element_3d_length(&bar3) - 13.0) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_truss_elements - 3D\n");
    }

    // End
    return rst;
}

bool c_test_discrete_elements()
{
    // Local Variables
    const double tol = 1.0e-12;
    bool rst;
    int i;
    double m[36], c[36], k[36], ke[16], mans[36], cans[36], kans[36];
    double u[4] = {0.0, 0.0, 0.0, 0.25};
    c_node nodes[3];
    c_spring_element_2d springs[2];
    c_damper_element_2d dampers[2];
    c_mass_element_2d masses[2];
    c_spring_element_3d spring3;
    c_node nodes3[2];
    double k3[36];

    // Initialization
    rst = true;
    nodes[0] = make_node(1, 2, 0.0, 0.0, 0.0);
    nodes[1] = make_node(2, 2, 0.0, 1.0, 0.0);
    nodes[2] = make_node(3, 2, 0.0, 2.0, 0.0);
    for (i = 0; i < 2; ++i)
    {
        springs[i].stiffness = (i == 0) ? 100.0 : 200.0;
        springs[i].node_1 = nodes[i];
        springs[i].node_2 = nodes[i + 1];
        springs[i].use_direction = false;
        dampers[i].damping_coefficient = (i == 0) ? 3.0 : 5.0;
        dampers[i].node_1 = nodes[i];
        dampers[i].node_2 = nodes[i + 1];
        dampers[i].use_direction = false;
        masses[i].mass = (i == 0) ? 1.5 : 2.5;
        masses[i].node_1 = nodes[i + 1];
    }

    // Expected results; only the y-direction DOFs (1, 3, 5) are coupled
    zero_matrix(6, 6, mans, 6);
    zero_matrix(6, 6, cans, 6);
    zero_matrix(6, 6, kans, 6);
    mans[INDEX(2, 2, 6)] = 1.5;
    mans[INDEX(3, 3, 6)] = 1.5;
    mans[INDEX(4, 4, 6)] = 2.5;
    mans[INDEX(5, 5, 6)] = 2.5;
    kans[INDEX(1, 1, 6)] = 100.0;
    kans[INDEX(1, 3, 6)] = -100.0;
    kans[INDEX(3, 1, 6)] = -100.0;
    kans[INDEX(3, 3, 6)] = 300.0;
    kans[INDEX(3, 5, 6)] = -200.0;
    kans[INDEX(5, 3, 6)] = -200.0;
    kans[INDEX(5, 5, 6)] = 200.0;
    cans[INDEX(1, 1, 6)] = 3.0;
    cans[INDEX(1, 3, 6)] = -3.0;
    cans[INDEX(3, 1, 6)] = -3.0;
    cans[INDEX(3, 3, 6)] = 8.0;
    cans[INDEX(3, 5, 6)] = -5.0;
    cans[INDEX(5, 3, 6)] = -5.0;
    cans[INDEX(5, 5, 6)] = 5.0;

    // Assemble the complete system
    c_assemble_discrete_system_2d(6, 2, masses, 2, dampers, 2, springs, 3,
        nodes, m, 6, c, 6, k, 6);
    if (!compare_matrices(6, 6, m, 6, mans, 6, tol) ||
        !compare_matrices(6, 6, c, 6, cans, 6, tol) ||
        !compare_matrices(6, 6, k, 6, kans, 6, tol) ||
        !c_is_symmetric(6, 6, k, 6))
    {
        rst = false;
        printf("TEST FAILED: c_test_discrete_elements - assembly\n");
    }

    // Assemble without dampers
    c_assemble_discrete_system_2d(6, 2, masses, 0, NULL, 2, springs, 3,
        nodes, m, 6, c, 6, k, 6);
    zero_matrix(6, 6, cans, 6);
    if (!compare_matrices(6, 6, c, 6, cans, 6, 0.0) ||
        !compare_matrices(6, 6, k, 6, kans, 6, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_discrete_elements - no dampers\n");
    }

    // Individual element routines
    c_spring_element_2d_stiffness_matrix(&springs[0], ke, 4);
    c_mass_element_2d_mass_matrix(&masses[0], m, 2);
    if (fabs(ke[INDEX(1, 1, 4)] - 100.0) > tol ||
        fabs(ke[INDEX(1, 3, 4)] + 100.0) > tol ||
        fabs(c_spring_element_2d_force(&springs[0], u) - 25.0) > tol ||
        fabs(c_damper_element_2d_force(&dampers[0], u) - 0.75) > tol ||
        fabs(m[INDEX(0, 0, 2)] - 1.5) > tol ||
        fabs(m[INDEX(1, 1, 2)] - 1.5) > tol ||
        fabs(m[INDEX(0, 1, 2)]) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_discrete_elements - elements\n");
    }

    // 3D zero-length spring with a user-defined axis
    nodes3[0] = make_node(1, 3, 0.0, 0.0, 0.0);
    nodes3[1] = make_node(2, 3, 0.0, 0.0, 0.0);
    spring3.stiffness = 14.0;
    spring3.node_1 = nodes3[0];
    spring3.node_2 = nodes3[1];
    spring3.direction[0] = 0.0;
    spring3.direction[1] = 0.0;
    spring3.direction[2] = 3.0;
    spring3.use_direction = true;
    c_spring_element_3d_stiffness_matrix(&spring3, k3, 6);
    if (fabs(k3[INDEX(2, 2, 6)] - 14.0) > tol ||
        fabs(k3[INDEX(2, 5, 6)] + 14.0) > tol ||
        fabs(k3[INDEX(0, 0, 6)]) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_discrete_elements - 3D user axis\n");
    }

    // End
    return rst;
}

bool c_test_generalized_alpha_integrator()
{
    // Local Variables
    const int nsteps = 200;
    const double dt = 1.0e-2;
    const double tol = 1.0e-3;
    bool rst;
    int i;
    double m = 1.0, c = 0.0, k = 4.0, f = 0.0;
    double u, v, a, us, vs, as, t;
    double *forces;
    c_structural_integrator obj;

    // Initialization
    rst = true;
    obj = c_create_dense_generalized_alpha_integrator(1, &m, 1, &c, 1, &k, 1,
        1.0);
    if (!obj) return false;

    // Step-by-step free vibration: u(t) = cos(2t)
    u = 1.0;
    v = 0.0;
    a = -4.0;
    for (i = 0; i < nsteps; ++i)
    {
        c_structural_integrator_step(obj, 1, &f, &f, dt, &u, &v, &a);
    }
    t = nsteps * dt;
    if (fabs(u - cos(2.0 * t)) > tol || fabs(v + 2.0 * sin(2.0 * t)) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_generalized_alpha_integrator - step\n");
    }

    // Force-history solve should reproduce the stepped result
    forces = (double*)calloc((size_t)(nsteps + 1), sizeof(double));
    if (!forces)
    {
        c_free_structural_integrator(obj);
        return false;
    }
    us = 1.0;
    vs = 0.0;
    as = -4.0;
    c_structural_integrator_solve(obj, 1, nsteps + 1, forces, 1, dt, &us, &vs,
        &as);
    if (fabs(us - u) > 1.0e-12 || fabs(vs - v) > 1.0e-12 ||
        fabs(as - a) > 1.0e-12)
    {
        rst = false;
        printf("TEST FAILED: c_test_generalized_alpha_integrator - solve\n");
    }

    // End
    free(forces);
    c_free_structural_integrator(obj);
    return rst;
}

static void oscillator(int n, double t, const double *x, double *dxdt,
    void *user_data)
{
    const double wn = *(const double*)user_data;
    dxdt[0] = x[1];
    dxdt[1] = -wn * wn * x[0];
    dxdt[2] = 1.0;
}

static void identity_coordinates(int n, double t, const double *x,
    double coordinates[3], void *user_data)
{
    coordinates[0] = x[0];
    coordinates[1] = x[1];
    coordinates[2] = t;
}

bool c_test_poincare_map_ode()
{
    // Local Variables
    const double pi = 3.14159265358979323846;
    const double tol = 1.0e-4;
    const int nbuffer = 10;
    bool rst;
    int i, pass, nactual;
    double tspan[2] = {0.0, 10.0};
    double iv[3] = {1.0, 0.0, 0.0};
    double pt[3] = {0.0, 0.0, 0.0};
    double nrm[3] = {0.0, 1.0, 0.0};
    double xb[10], yb[10], zb[10];
    double wn = 1.0;
    c_plane pln;

    // Initialization
    rst = true;
    c_plane_from_point_and_normal(pt, nrm, &pln);

    // Crossings of v = 0 occur at t = pi, 2 pi, and 3 pi
    for (pass = 0; pass < 2; ++pass)
    {
        c_poincare_map_ode(oscillator, tspan, 3, iv, 10001, &pln,
            DYN_POINCARE_TWO_SIDED, DYN_RUNGE_KUTTA_45,
            (pass == 0) ? NULL : identity_coordinates, nbuffer, xb, yb, zb,
            &nactual, &wn);
        if (nactual != 3)
        {
            rst = false;
            printf("TEST FAILED: c_test_poincare_map_ode - count (%d)\n",
                pass);
            continue;
        }
        for (i = 0; i < 3; ++i)
        {
            if (fabs(xb[i] - ((i % 2 == 0) ? -1.0 : 1.0)) > tol ||
                fabs(yb[i]) > tol ||
                fabs(zb[i] - (i + 1) * pi) > tol)
            {
                rst = false;
                printf("TEST FAILED: c_test_poincare_map_ode - point (%d)\n",
                    pass);
                break;
            }
        }
    }

    // End
    return rst;
}

static void count_errors(int code, const char *message, void *user_data)
{
    int *count = (int*)user_data;
    if (message && message[0] != '\0') *count += 1;
}

bool c_test_error_reporting()
{
    // Local Variables
    bool rst;
    int count;
    char msg[16];
    double k[16], f = 0.0, u = 0.0, v = 0.0, a = 0.0;
    c_truss_element_2d bar;

    // Initialization
    rst = true;
    count = 0;
    bar.material.modulus = 1.0;
    bar.material.density = 1.0;
    bar.material.poissons_ratio = 0.3;
    bar.area = 1.0;
    bar.node_1 = make_node(1, 2, 0.0, 0.0, 0.0);
    bar.node_2 = make_node(2, 2, 1.0, 0.0, 0.0);
    c_clear_error();
    c_set_error_handler(count_errors, &count);

    // An undersized leading dimension must be reported, not terminate
    c_truss_element_2d_stiffness_matrix(&bar, k, 2);
    if (c_get_last_error() != DYN_INVALID_INPUT_ERROR || count != 1 ||
        c_get_last_error_message((int)sizeof(msg), msg) < (int)sizeof(msg) ||
        msg[sizeof(msg) - 1] != '\0')
    {
        rst = false;
        printf("TEST FAILED: c_test_error_reporting - leading dimension\n");
    }

    // A NULL handle must be reported
    c_structural_integrator_step(NULL, 1, &f, &f, 0.1, &u, &v, &a);
    if (c_get_last_error() != DYN_NULL_POINTER_ERROR || count != 2)
    {
        rst = false;
        printf("TEST FAILED: c_test_error_reporting - null handle\n");
    }

    // Clearing resets the error; valid calls do not record one
    c_clear_error();
    c_truss_element_2d_stiffness_matrix(&bar, k, 4);
    if (c_get_last_error() != DYN_NO_ERROR || count != 2 ||
        c_get_last_error_message((int)sizeof(msg), msg) != 0)
    {
        rst = false;
        printf("TEST FAILED: c_test_error_reporting - clear\n");
    }

    // End
    c_set_error_handler(NULL, NULL);
    return rst;
}

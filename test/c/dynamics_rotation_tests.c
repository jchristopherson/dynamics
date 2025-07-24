#include "dynamics_rotation_tests.h"
#include "dynamics_c_test_helper.h"
#include "dynamics.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

bool c_test_rotation_x()
{
    // Local Variables
    bool rst;
    double angle, r[9], ans[9];
    const double tol = 1.0e-8;
    
    // Initialization
    rst = true;
    angle = ((double)rand()) / RAND_MAX;
    ans[0] = 1.0;
    ans[1] = 0.0;
    ans[2] = 0.0;
    ans[3] = 0.0;
    ans[4] = cos(angle);
    ans[5] = sin(angle);
    ans[6] = 0.0;
    ans[7] = -sin(angle);
    ans[8] = cos(angle);

    // Test
    c_rotate_x(angle, r, 3);
    if (!compare_matrices(3, 3, r, 3, ans, 3, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_rotation_x -1\n");
    }

    // End
    return rst;
}


bool c_test_rotation_y()
{
    // Local Variables
    bool rst;
    double angle, r[9], ans[9];
    const double tol = 1.0e-8;
    
    // Initialization
    rst = true;
    angle = ((double)rand()) / RAND_MAX;
    ans[0] = cos(angle);
    ans[1] = 0.0;
    ans[2] = -sin(angle);
    ans[3] = 0.0;
    ans[4] = 1.0;
    ans[5] = 0.0;
    ans[6] = sin(angle);
    ans[7] = 0.0;
    ans[8] = cos(angle);

    // Test
    c_rotate_y(angle, r, 3);
    if (!compare_matrices(3, 3, r, 3, ans, 3, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_rotation_y -1\n");
    }

    // End
    return rst;
}


bool c_test_rotation_z()
{
    // Local Variables
    bool rst;
    double angle, r[9], ans[9];
    const double tol = 1.0e-8;
    
    // Initialization
    rst = true;
    angle = ((double)rand()) / RAND_MAX;
    ans[0] = cos(angle);
    ans[1] = sin(angle);
    ans[2] = 0.0;
    ans[3] = -sin(angle);
    ans[4] = cos(angle);
    ans[5] = 0.0;
    ans[6] = 0.0;
    ans[7] = 0.0;
    ans[8] = 1.0;

    // Test
    c_rotate_z(angle, r, 3);
    if (!compare_matrices(3, 3, r, 3, ans, 3, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_rotation_z -1\n");
    }

    // End
    return rst;
}


bool c_test_rotation()
{
    // Local Variables
    bool rst;
    double tx, rx[9], r[9], i[3], j[3], k[3], ix[3], jx[3], kx[3];
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    tx = ((double)rand()) / RAND_MAX;
    c_rotate_x(tx, rx, 3);
    i[0] = 1.0;     i[1] = 0.0;     i[2] = 0.0;
    j[0] = 0.0;     j[1] = 1.0;     j[2] = 0.0;
    k[0] = 0.0;     k[1] = 0.0;     k[2] = 1.0;
    mtx_mult(3, 1, 3, rx, 3, i, 3, ix, 3);
    mtx_mult(3, 1, 3, rx, 3, j, 3, jx, 3);
    mtx_mult(3, 1, 3, rx, 3, k, 3, kx, 3);

    // Test
    c_rotate(ix, jx, kx, r, 3);
    if (!compare_matrices(3, 3, rx, 3, r, 3, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_rotation -1\n");
    }

    // End
    return rst;
}


bool c_test_acceleration_transform()
{
    // Local Variables
    bool rst;
    double alphavec[3], omegavec[3], alpha[9], w[9], wt[9], a[3], x[3], A[16], 
        ans[16], wwt[9], awwtx[9];
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(3, alphavec);
    create_random_vector(3, omegavec);
    to_skew_symmetric_matrix(alphavec, alpha, 3);
    to_skew_symmetric_matrix(omegavec, w, 3);
    transpose(3, 3, w, 3, wt, 3);
    zero_matrix(4, 4, ans, 4);

    // Compute w * w^T
    mtx_mult(3, 3, 3, w, 3, wt, 3, wwt, 3);

    // Compute alpha - w * w^T
    mtx_subtract(3, 3, alpha, 3, wwt, 3, ans, 4);

    // Compute (alpha - w * w^T) * x
    mtx_mult(3, 1, 3, ans, 4, x, 3, awwtx, 3);

    // Compute a - (alpha - w * w^T) * x
    vec_subtract(3, a, awwtx, &ans[12]);

    // Test
    c_acceleration_transform(alphavec, omegavec, a, x, A, 4);
    if (!compare_matrices(4, 4, ans, 4, A, 4, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_acceleration_transform -1\n");
    }

    // End
    return rst;
}


bool c_test_velocity_transform()
{
    // Local Variables
    bool rst;
    double omega[3], v[3], x[3], V[16], ans[16], wx[3];
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(3, omega);
    create_random_vector(3, v);
    create_random_vector(3, x);
    zero_matrix(4, 4, ans, 4);
    to_skew_symmetric_matrix(omega, ans, 4);

    // Compute omega * x
    mtx_mult(3, 1, 3, ans, 4, x, 3, wx, 3);

    // Compute v - omega * x
    mtx_subtract(3, 1, v, 3, wx, 3, &ans[12], 3);

    // Test
    c_velocity_transform(omega, v, x, V, 4);
    if (!compare_matrices(4, 4, V, 4, ans, 4, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_velocity_transform -1\n");
    }

    // End
    return rst;
}


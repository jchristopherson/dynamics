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

bool c_test_quaternion_from_array()
{
    // Local Variables
    bool rst;
    double x[4];
    c_quaternion q;
    const double tol = 1.0e-8;

    create_random_vector(4, x);
    c_quaternion_from_array(x, &q);
    rst = compare_quaternion_to_array(x, &q, tol);
    if (!rst)
    {
        printf("TEST FAILED: c_test_quaternion_from_array -1\n");
    }
    return rst;
}

bool c_test_quaternion_from_angle_axis()
{
    // Local Variables
    bool rst;
    double angle, axis[3], ans[4];
    c_quaternion q;
    const double tol = 1.0e-8;

    // Initialization
    create_random_vector(3, axis);
    angle = 0.75;
    ans[0] = cos(0.5 * angle);
    ans[1] = sin(0.5 * angle) * axis[0];
    ans[2] = sin(0.5 * angle) * axis[1];
    ans[3] = sin(0.5 * angle) * axis[2];
    c_quaternion_from_angle_axis(angle, axis, &q);

    // Test
    rst = compare_quaternion_to_array(ans, &q, tol);
    if (!rst)
    {
        printf("TEST FAILED: c_test_quaternion_from_angle_axis -1\n");
    }
    return rst;
}

bool c_test_quaternion_from_matrix()
{
    // Local Variables
    bool rst;
    double Rx[9], ax, axis[3];
    c_quaternion q, p;
    const double tol = 1.0e-8;
    
    // Initialization
    ax = 0.25;
    axis[0] = 1.0;
    axis[1] = 0.0;
    axis[2] = 0.0;
    c_rotate_x(ax, Rx, 3);
    c_quaternion_from_matrix(Rx, 3, &q);
    c_quaternion_from_angle_axis(ax, axis, &p);

    c_quaternion_normalize(&q);
    c_quaternion_normalize(&p);

    // Test
    rst = compare_quaternions(&q, &p, tol);
    if (!rst)
    {
        printf("TEST FAILED: c_test_quaternion_from_matrix -1\n");
    }

    // End
    return rst;
}

bool c_test_quaternion_normalize()
{
    // Local Variables
    bool rst;
    c_quaternion q;
    double x[4], xnorm;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(4, x);
    c_quaternion_from_array(x, &q);
    xnorm = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
    x[0] /= xnorm;
    x[1] /= xnorm;
    x[2] /= xnorm;
    x[3] /= xnorm;
    c_quaternion_normalize(&q);

    // Test
    rst = compare_quaternion_to_array(x, &q, tol);
    if (!rst)
    {
        printf("TEST FAILED: test_quaternion_normalize -1\n");
    }
    return rst;
}

bool c_test_quaternion_add()
{
    // Local Variables
    bool rst;
    c_quaternion q1, q2, q;
    double x1[4], x2[4], x[4];
    int i;
    const double tol = 1.0e-8;

    // Initialization
    create_random_vector(4, x1);
    create_random_vector(4, x2);
    for (i = 0; i < 4; ++i) x[i] = x1[i] + x2[i];
    c_quaternion_from_array(x1, &q1);
    c_quaternion_from_array(x2, &q2);
    c_quaternion_add(&q1, &q2, &q);

    // Test
    rst = compare_quaternion_to_array(x, &q, tol);
    if (!rst)
    {
        printf("TEST FAILED: c_test_quaternion_add -1\n");
    }
    return rst;
}

bool c_test_quaternion_subtract()
{
    // Local Variables
    bool rst;
    c_quaternion q1, q2, q;
    double x1[4], x2[4], x[4];
    int i;
    const double tol = 1.0e-8;

    // Initialization
    create_random_vector(4, x1);
    create_random_vector(4, x2);
    for (i = 0; i < 4; ++i) x[i] = x1[i] - x2[i];
    c_quaternion_from_array(x1, &q1);
    c_quaternion_from_array(x2, &q2);
    c_quaternion_subtract(&q1, &q2, &q);

    // Test
    rst = compare_quaternion_to_array(x, &q, tol);
    if (!rst)
    {
        printf("TEST FAILED: c_test_quaternion_subtract -1\n");
    }
    return rst;
}

bool c_test_quaternion_multiply()
{
    // Local Variables
    bool rst1, rst2, rst;
    c_quaternion q, qc, qrb, qq;
    double R[9], rb[3], rg[3], rgq[3], angle, axis[3], rbb[4], rgg[3];
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    angle = -0.75;
    axis[0] = 1.0;  axis[1] = 0.0;  axis[2] = 0.0;
    c_rotate_x(angle, R, 3);
    create_random_vector(3, rb);
    rbb[0] = 0.0; rbb[1] = rb[0]; rbb[2] = rb[1]; rbb[3] = rb[2];
    c_quaternion_from_angle_axis(angle, axis, &q);
    c_quaternion_from_array(rbb, &qrb);
    mtx_mult(3, 1, 3, R, 3, rb, 3, rg, 3);  // rg = R * rb
    c_quaternion_normalize(&q);
    c_quaternion_conjugate(&q, &qc);

    // Compute rg = q * rb * q*
    c_quaternion_multiply(&q, &qrb, &qq);
    c_quaternion_multiply(&qq, &qc, &qrb);  // overwrite qrb
    rgq[0] = qrb.x; rgq[1] = qrb.y; rgq[2] = qrb.z;

    // Test
    rst1 = compare_arrays(3, rgq, rg, tol);
    if (!rst1)
    {
        printf("TEST FAILED: c_test_quaternion_multiply -1\n");
        rst = false;
    }

    // Test 2 - c_quaternion_rotate
    c_quaternion_rotate(&q, rb, rgg);
    rst2 = compare_arrays(3, rg, rgg, tol);
    if (!rst2)
    {
        printf("TEST FAILED: c_test_quaternion_multiply -2\n");
        rst = false;
    }

    return rst;
}

bool c_test_quaternion_divide()
{
    // Local Variables
    bool rst1, rst2, rst;
    c_quaternion q1, q2, q2c, q, qi, qa, q2inv;
    double x1[4], x2[4], q2abs;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(4, x1);
    create_random_vector(4, x2);
    c_quaternion_from_array(x1, &q1);
    c_quaternion_from_array(x2, &q2);
    c_quaternion_conjugate(&q2, &q2c);
    q2abs = c_quaternion_abs(&q2);
    c_quaternion_scale(1.0 / (q2abs * q2abs), &q2c, &qi); // q2* / |q2|^2
    c_quaternion_multiply(&q1, &qi, &qa);    // q1 * q2* / |q2|^2
    c_quaternion_divide(&q1, &q2, &q);

    // Test 1
    rst1 = compare_quaternions(&q, &qa, tol);
    if (!rst1)
    {
        printf("TEST FAILED: c_test_quaternion_divide -1\n");
        rst = false;
    }

    // Test 2
    c_quaternion_inverse(&q2, &q2inv);
    rst2 = compare_quaternions(&qi, &q2inv, tol);
    if (!rst2)
    {
        printf("TEST FAILED: c_test_quaternion_divide -2\n");
        rst = false;
    }

    // End
    return rst;
}

bool c_test_quaternion_abs()
{
    // Local Variables
    bool rst;
    c_quaternion q;
    double x[4], xabs, qabs;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    create_random_vector(4, x);
    c_quaternion_from_array(x, &q);
    xabs = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
    qabs = c_quaternion_abs(&q);

    // Test
    if (fabs(xabs - qabs) > tol)
    {
        printf("TEST FAILED: c_test_quaternion_abs -1\n");
        rst = false;
    }
    return rst;
}

bool c_test_quaternion_to_matrix()
{
    // Local Variables
    bool rst;
    c_quaternion q;
    double R[9], angle, axis[3], Rq[9];
    const double tol = 1.0e-8;

    // Initialization
    angle = -1.25;
    axis[0] = 1.0;  axis[1] = 0.0;  axis[2] = 0.0;
    c_quaternion_from_angle_axis(angle, axis, &q);
    c_rotate_x(angle, R, 3);
    c_quaternion_to_matrix(&q, Rq, 3);

    // Test
    rst = compare_matrices(3, 3, R, 3, Rq, 3, tol);
    if (!rst)
    {
        printf("TEST FAILED: c_test_quaternion_to_matrix -1\n");
    }
    return rst;
}

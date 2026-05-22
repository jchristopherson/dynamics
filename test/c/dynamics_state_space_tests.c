#include <math.h>
#include <complex.h>
#include <stdio.h>
#include "dynamics.h"
#include "dynamics_c_test_helper.h"
#include "dynamics_state_space_tests.h"

bool c_test_state_space_initialize()
{
    // Local Variables
    const double tol = 1.0e-8;
    const double m1 = 1.5;
    const double m2 = 3.2;
    const double b1 = 0.35;
    const double b2 = 4.25;
    const double k1 = 125.5;
    const double k2 = 250.1;
    const int n_out = 1;
    const int n_in = 2;
    bool rst;
    double M[4], B[4], K[4], Aans[16], Bans[8], Cans[4], Dans[2];
    c_state_space_model mdl;

    // Initialization
    rst = true;
    M[0] = m1;      M[1] = 0.0;     M[2] = 0.0;     M[3] = m2;
    B[0] = b1 + b2; B[1] = -b2;     B[2] = -b2;     B[3] = b2;
    K[0] = k1 + k2; K[1] = -k2;     K[2] = -k2;     K[3] = k2;

    Aans[0] = 0.0;   Aans[1] = 0.0;   Aans[2] = -K[0]/m1;  Aans[3] = k2/m2;
    Aans[4] = 0.0;   Aans[5] = 0.0;   Aans[6] = k2/m1;     Aans[7] = -k2/m2;
    Aans[8] = 1.0;   Aans[9] = 0.0;   Aans[10] = -B[0]/m1; Aans[11] = b2/m2;
    Aans[12] = 0.0;  Aans[13] = 1.0;  Aans[14] = b2/m1;    Aans[15] = -b2/m2;

    Bans[0] = 0.0;  Bans[1] = 0.0;  Bans[2] = 1.0/m1;   Bans[3] = 0.0;
    Bans[4] = 0.0;  Bans[5] = 0.0;  Bans[6] = 0.0;      Bans[7] = 1.0 / m2;

    Cans[0] = 1.0;  Cans[1] = 1.0;  Cans[2] = 1.0;  Cans[3] = 1.0;
    Dans[0] = 0.0;  Dans[1] = 0.0;

    // Test
    c_create_state_space_model(4, 1, M, 2, B, 2, K, 2, &mdl);
    if (mdl.dimension != 4)
    {
        rst = false;
        printf("TEST FAILED: c_test_state_space_initialize -1\n");
    }
    if (mdl.n_inputs != 2)
    {
        rst = false;
        printf("TEST FAILED: c_test_state_space_initialize -2\n");
    }
    if (mdl.n_outputs != 1)
    {
        rst = false;
        printf("TEST FAILED: c_test_state_space_initialize -3\n");
    }
    if (!rst) 
    {
        c_free_state_space_model(&mdl);
        return false;
    }

    if (!compare_matrices(4, 4, mdl.A, 4, Aans, 4, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_state_space_initialize -4\n");
    }
    if (!compare_matrices(4, 2, mdl.B, 4, Bans, 4, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_state_space_initialize -5\n");
    }
    if (!compare_matrices(1, 4, mdl.C, 1, Cans, 1, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_state_space_initialize -6\n");
    }
    if (!compare_matrices(1, 2, mdl.D, 1, Dans, 1, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_state_space_initialize -7\n");
    }

    // End
    c_free_state_space_model(&mdl);
    return rst;
}

bool c_test_state_space_poles_zeros()
{
    const double tol = 1.0e-8;
    const double m = 0.5;
    const double b = 1.65;
    const double k = 1.0e5;
    bool rst;
    c_transfer_function tf;
    c_state_space_model mdl;
    double complex z[1], p[2], zans[1], pans[2];
    int n;

    // Initialization
    rst = true;
    c_alloc_transfer_function(1, 2, &tf);
    tf.numerator.coefficients[0] = k;
    tf.numerator.coefficients[1] = b;
    tf.denominator.coefficients[0] = k;
    tf.denominator.coefficients[1] = b;
    tf.denominator.coefficients[2] = m;
    c_to_ccf_state_space(&tf, &mdl);
    c_transfer_function_poles(&tf, 2, pans);
    c_transfer_function_zeros(&tf, 1, zans);

    // Test
    c_state_space_poles(&mdl, 2, p);
    if (!compare_complex_arrays(2, p, pans, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_state_space_poles_zeros -1\n");
    }

    c_state_space_zeros(&mdl, 1, z, &n);
    if (n != 1)
    {
        rst = false;
        printf("TEST FAILED: c_test_state_space_poles_zeros -2\n");
        c_free_transfer_function(&tf);
        c_free_state_space_model(&mdl);
        return rst;
    }
    if (!compare_complex_arrays(1, z, zans, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_state_space_poles_zeros -3\n");
    }

    // End
    c_free_transfer_function(&tf);
    c_free_state_space_model(&mdl);
    return rst;
}

bool c_test_state_space_to_transfer_function()
{
    const double tol = 1.0e-8;
    const double m = 0.5;
    const double b = 1.65;
    const double k = 1.0e5;
    const double omega = 500.0;
    bool rst;
    double complex s, z, zans;
    c_transfer_function tf;
    c_state_space_model mdl;

    // Initialization
    rst = true;
    c_alloc_transfer_function(1, 2, &tf);
    tf.numerator.coefficients[0] = k;
    tf.numerator.coefficients[1] = b;
    tf.denominator.coefficients[0] = k;
    tf.denominator.coefficients[1] = b;
    tf.denominator.coefficients[2] = m;
    c_to_ccf_state_space(&tf, &mdl);
    s = omega * I;
    c_evaluate_transfer_function(&tf, 1, &s, &zans);

    // Test
    c_state_space_transfer_function(&mdl, mdl.n_inputs, mdl.n_outputs, 1, &s, 
        &z, 1);
    if (!compare_complex_arrays(1, &z, &zans, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_state_space_to_transfer_function -1\n");
    }

    // End
    c_free_transfer_function(&tf);
    c_free_state_space_model(&mdl);
    return rst;
}

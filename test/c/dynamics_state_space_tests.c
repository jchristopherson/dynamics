#include <math.h>
#include <complex.h>
#include <stdio.h>
#include "dynamics.h"
#include "dynamics_c_test_helper.h"
#include "dynamics_state_space_tests.h"

#define DISCRETE_CHECK(condition) \
    do { \
        if (!(condition)) { \
            printf("Discrete C test failed at line %d: %s\n", __LINE__, #condition); \
            goto cleanup; \
        } \
    } while (0)

static bool test_discrete_models(void)
{
    const double ts = 0.1, tol = 1e-10;
    double a = -2, b = 1, c = 3, d = 0.5, u = 2, x = 4, next, y;
    double delay_coeff[] = {0, 0, 2}, gain = 2, one = 1;
    double omega = acos(-1) / (2 * ts);
    double complex roots[1], response;
    int count;
    bool ok = false;
    c_state_space_model plant = {1, 1, 1, &a, &b, &c, &d};
    c_state_space_model disc = {0}, restored = {0}, realization = {0};
    c_transfer_function tf = {0};
    c_transfer_function delay = {{2, delay_coeff}, {0, &one}};
    c_transfer_function constant = {{0, &gain}, {0, &one}};
    c_clear_error();

    for (int method = DYN_DISCRETIZE_ZOH; method <= DYN_DISCRETIZE_TUSTIN; ++method) {
        DISCRETE_CHECK(c_c2d(&plant, ts, method, &disc) == 0);
        DISCRETE_CHECK(c_d2c(&disc, ts, method, &restored) == 0);
        DISCRETE_CHECK(fabs(restored.A[0] - a) < tol && fabs(restored.B[0] - b) < tol &&
            fabs(restored.C[0] - c) < tol && fabs(restored.D[0] - d) < tol);
        c_free_state_space_model(&restored);
        if (method == DYN_DISCRETIZE_ZOH) {
            const double pole = exp(a * ts), bd = (pole - 1) * b / a;
            DISCRETE_CHECK(c_discrete_state_space_step(&disc, ts, &u, &x, &next, &y) == 0);
            DISCRETE_CHECK(fabs(next - (pole * x + bd * u)) < tol && fabs(y - 13) < tol);
            DISCRETE_CHECK(c_discrete_state_space_roots(&disc, ts, false, 1, roots, &count) == 0);
            DISCRETE_CHECK(count == 1 && cabs(roots[0] - pole) < tol);
            DISCRETE_CHECK(c_discrete_state_space_roots(&disc, ts, true, 0, NULL, &count) == 0);
            DISCRETE_CHECK(count == 1);
            DISCRETE_CHECK(c_discrete_state_space_roots(&disc, ts, true, 0, roots, &count) == -1);
            DISCRETE_CHECK(count == 1 && c_get_last_error() == DYN_ARRAY_SIZE_ERROR);
            DISCRETE_CHECK(c_discrete_state_space_roots(&disc, ts, true, 1, roots, &count) == 0);
            DISCRETE_CHECK(cabs(roots[0] - (pole - c * bd / d)) < tol);
            DISCRETE_CHECK(c_discrete_state_space_to_transfer_function(&disc, ts, &tf) == 0);
            DISCRETE_CHECK(c_discrete_transfer_function_roots(&tf, ts, true, 1, roots, &count) == 0);
            DISCRETE_CHECK(count == 1 && cabs(roots[0] - (pole - c * bd / d)) < tol);
        }
        c_free_state_space_model(&disc);
    }

    DISCRETE_CHECK(c_evaluate_discrete_transfer_function(&delay, ts, 1, &omega, &response) == 0);
    DISCRETE_CHECK(cabs(response + 2) < tol);
    for (int observable = 0; observable < 2; ++observable) {
        double input = 1, state[2] = {0}, updated[2];
        DISCRETE_CHECK(c_discrete_transfer_function_to_state_space(&delay, ts, observable, &realization) == 0);
        DISCRETE_CHECK(realization.dimension == 2);
        for (int sample = 0; sample < 3; ++sample) {
            DISCRETE_CHECK(c_discrete_state_space_step(&realization, ts, &input, state, updated, &y) == 0);
            DISCRETE_CHECK(fabs(y - (sample == 2 ? 2 : 0)) < tol);
            state[0] = updated[0];
            state[1] = updated[1];
            input = 0;
        }
        c_free_state_space_model(&realization);
    }
    DISCRETE_CHECK(c_discrete_transfer_function_to_state_space(&constant, ts, false, &realization) == 0);
    DISCRETE_CHECK(realization.dimension == 0);
    DISCRETE_CHECK(c_discrete_state_space_step(&realization, ts, &u, NULL, NULL, &y) == 0);
    DISCRETE_CHECK(fabs(y - 4) < tol);
    c_free_transfer_function(&tf);
    DISCRETE_CHECK(c_discrete_state_space_to_transfer_function(&realization, ts, &tf) == 0);
    DISCRETE_CHECK(tf.numerator.order == 0 && fabs(tf.numerator.coefficients[0] - gain) < tol);

    DISCRETE_CHECK(c_c2d(&plant, NAN, DYN_DISCRETIZE_ZOH, &disc) == -1);
    DISCRETE_CHECK(c_get_last_error() == DYN_INVALID_INPUT_ERROR);
    a = 2 / ts;
    DISCRETE_CHECK(c_c2d(&plant, ts, DYN_DISCRETIZE_TUSTIN, &disc) == -1);
    a = -1;
    DISCRETE_CHECK(c_d2c(&plant, ts, DYN_DISCRETIZE_TUSTIN, &restored) == -1);
    DISCRETE_CHECK(c_d2c(&plant, ts, DYN_DISCRETIZE_ZOH, &restored) == -1);
    DISCRETE_CHECK(c_get_last_error() == DYN_INVALID_INPUT_ERROR);
    ok = true;
cleanup:
    c_free_transfer_function(&tf);
    c_free_state_space_model(&disc);
    c_free_state_space_model(&restored);
    c_free_state_space_model(&realization);
    c_clear_error();
    return ok;
}

static bool test_discrete_feedback(void)
{
    const double ts = 0.1, tol = 1e-10;
    double a[] = {0.9, 0, 0.1, 0.8}, b[] = {0, 0.2}, c[] = {1, 0}, d = 0.1;
    double k[] = {1, 0.5}, ki = 0.3, bad_k[] = {0, 0}, bad_ki = 0;
    double x[] = {1, 2, 3}, next[3], y[2], reference = 2;
    double ise0, ise1, cost;
    bool stable, ok = false;
    c_discrete_optimizer_options options = {1e-9, 250, 0.01};
    c_state_space_model plant = {2, 1, 1, a, b, c, &d}, loop = {0};

    DISCRETE_CHECK(c_discrete_integral_feedback_apply(&plant, ts, k, &ki, &loop) == 0);
    DISCRETE_CHECK(loop.dimension == 3 && loop.n_inputs == 1 && loop.n_outputs == 2);
    DISCRETE_CHECK(c_discrete_state_space_step(&loop, ts, &reference, x, next, y) == 0);
    DISCRETE_CHECK(fabs(y[0] + 2.9) < tol && fabs(y[1] - 0.71) < tol);
    DISCRETE_CHECK(fabs(next[0] - 1.1) < tol && fabs(next[1] - 1.02) < tol);
    DISCRETE_CHECK(fabs(next[2] - (x[2] + ts * (y[1] - reference))) < tol);
    DISCRETE_CHECK(c_discrete_integral_feedback_analyze(&plant, ts, k, &ki, &stable, &ise0) == 0);
    DISCRETE_CHECK(stable && isfinite(ise0));
    DISCRETE_CHECK(c_discrete_integral_feedback_optimize(&plant, ts, k, &ki, NULL, &cost) == 0);
    DISCRETE_CHECK(c_discrete_integral_feedback_analyze(&plant, ts, k, &ki, &stable, &ise1) == 0);
    DISCRETE_CHECK(stable && ise1 < ise0 && fabs(cost - ise1) < 1e-6);
    DISCRETE_CHECK(c_discrete_integral_feedback_optimize(&plant, ts, k, &ki, &options, &cost) == 0);
    DISCRETE_CHECK(c_discrete_integral_feedback_analyze(&plant, ts, k, &ki, NULL, &ise1) == 0);
    DISCRETE_CHECK(isfinite(ise1) && cost > ise1);
    DISCRETE_CHECK(c_discrete_integral_feedback_optimize(&plant, ts, bad_k, &bad_ki, NULL, NULL) == -1);
    DISCRETE_CHECK(c_get_last_error() == DYN_INVALID_INPUT_ERROR);
    ok = true;
cleanup:
    c_free_state_space_model(&loop);
    c_clear_error();
    return ok;
}

bool c_test_discrete_controls(void)
{
    return test_discrete_models() && test_discrete_feedback();
}

#undef DISCRETE_CHECK

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
    c_create_state_space_model(2, 1, M, 2, B, 2, K, 2, &mdl);
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

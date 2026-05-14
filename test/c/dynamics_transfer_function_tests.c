#include "dynamics_c_test_helper.h"
#include "dynamics_transfer_function_tests.h"
#include <stdio.h>
#include <complex.h>
#include <math.h>

bool c_test_tf_evaluate()
{
    // Local Variables
    const double tol = 1.0e-8;
    const int n = 10;
    const double m = 1.5;
    const double b = 25.0;
    const double k = 125.3;
    bool rst;
    int i;
    double omega[n];
    double complex s[n], z[n], ans[n];
    c_transfer_function tf;

    // Initialization
    rst = true;
    create_random_vector(n, omega);
    for (i = 0; i < n; ++i) s[i] = omega[i] * I;
    c_alloc_transfer_function(0, 2, &tf);
    tf.numerator.coefficients[0] = 1.0;
    tf.denominator.coefficients[0] = k;
    tf.denominator.coefficients[1] = b;
    tf.denominator.coefficients[2] = m;

    // Compute the solution
    for (i = 0; i < n; ++i) ans[i] = 1.0 / (m * s[i] * s[i] + b * s[i] + k);

    // Test
    c_evaluate_transfer_function(&tf, n, s, z);
    if (!compare_complex_arrays(n, z, ans, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_tf_evaluate -1\n");
        printf("\nAnswer:\n");
        for (i = 0; i < n; ++i) 
            printf("%i\t%0.6f + %0.6fi\n", i + 1, creal(ans[i]), cimag(ans[i]));
        printf("\nComputed:\n");
        for (i = 0; i < n; ++i) 
            printf("%i.\t%0.6f + %0.6fi\n", i + 1, creal(z[i]), cimag(z[i]));
    }

    // End
    c_free_transfer_function(&tf);
    return rst;
}

bool c_test_ccf_form()
{
    // Local Variables
    const double tol = 1.0e-8;
    const int n = 4;
    const int m = 5;
    const int lda = n;
    const int ldb = n;
    const int ldc = 1;
    const int ldd = 1;
    bool rst;
    int i, j;
    double arg, x[m], y[n], A[n*n], B[n], C[n], D[1];
    c_transfer_function tf;
    c_state_space_model mdl;

    // Initialization
    rst = true;
    create_random_vector(m, x);
    create_random_vector(n, y);
    c_alloc_transfer_function(n - 1, m - 1, &tf);
    for (i = 0; i < n * n; ++i) A[i] = 0.0;
    for (i = 0; i < n; ++i) B[i] = 0.0;
    for (i = 0; i < n; ++i) C[i] = 0.0;
    D[0] = 0.0;
    for (i = 0; i < m; ++i) tf.denominator.coefficients[i] = x[i];
    for (i = 0; i < n; ++i) tf.numerator.coefficients[i] = y[i];

    // Define the solution
    arg = x[n];
    for (i = 0; i < n; ++i) y[i] = C[i] = y[i] / arg;
    for (i = 0; i < m; ++i) x[i] = x[i] / arg;
    
    for (i = 0; i < n; ++i) A[INDEX(n-1,i,lda)] = -x[i];
    B[INDEX(n-1,0,ldb)] = 1.0;
    for (i = 0; i < n - 1; ++i) A[INDEX(i,i+1,lda)] = 1.0;

    // Test
    c_to_ccf_state_space(&tf, &mdl);

    if (!compare_matrices(n, n, A, lda, mdl.A, n, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_ccf_form -1\n");
    }
    if (!compare_matrices(n, 1, B, ldb, mdl.B, n, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_ccf_form -2\n");
    }
    if (!compare_matrices(1, n, C, ldc, mdl.C, 1, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_ccf_form -3\n");
    }
    if (!compare_matrices(1, 1, D, ldd, mdl.D, 1, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_ccf_form -4\n");
    }

    // End
    c_free_transfer_function(&tf);
    c_free_state_space_model(&mdl);
    return rst;
}

bool c_test_ocf_form()
{
    // Local Variables
    const double tol = 1.0e-8;
    const int n = 4;
    const int m = 5;
    const int lda = n;
    const int ldb = n;
    const int ldc = 1;
    const int ldd = 1;
    bool rst;
    int i, j;
    double arg, x[m], y[n], A[n*n], B[n], C[n], D[1];
    c_transfer_function tf;
    c_state_space_model mdl;

    // Initialization
    rst = true;
    create_random_vector(m, x);
    create_random_vector(n, y);
    c_alloc_transfer_function(n - 1, m - 1, &tf);
    for (i = 0; i < n * n; ++i) A[i] = 0.0;
    for (i = 0; i < n; ++i) B[i] = 0.0;
    for (i = 0; i < n; ++i) C[i] = 0.0;
    D[0] = 0.0;
    for (i = 0; i < m; ++i) tf.denominator.coefficients[i] = x[i];
    for (i = 0; i < n; ++i) tf.numerator.coefficients[i] = y[i];

    // Define the solution
    arg = x[n];
    for (i = 0; i < n; ++i) y[i] = y[i] / arg;
    for (i = 0; i < m; ++i) x[i] = x[i] / arg;
    
    for (i = 0; i < n; ++i) A[INDEX(i,n-1,lda)] = -x[i];
    for (i = 0; i < n - 1; ++i) A[INDEX(i+1,i,lda)] = 1.0;
    for (i = 0; i < n; ++i) B[i] = y[i];
    C[INDEX(0,n-1,ldc)] = 1.0;

    // Test
    c_to_ocf_state_space(&tf, &mdl);

    if (!compare_matrices(n, n, A, lda, mdl.A, n, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_ocf_form -1\n");
    }
    if (!compare_matrices(n, 1, B, ldb, mdl.B, n, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_ocf_form -2\n");
    }
    if (!compare_matrices(1, n, C, ldc, mdl.C, 1, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_ocf_form -3\n");
    }
    if (!compare_matrices(1, 1, D, ldd, mdl.D, 1, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_ocf_form -4\n");
    }

    // End
    c_free_transfer_function(&tf);
    c_free_state_space_model(&mdl);
    return rst;
}

bool c_test_poles_zeros()
{
    // Local Variables
    const double tol = 1.0e-8;
    const double m = 1.5;
    const double b = 25.0;
    const double k = 125.3;
    const double complex p1 = (-b + sqrt(fabs(b*b - 4.0 * m * k))*I) / (2.0 * m);
    const double complex p2 = (-b - sqrt(fabs(b*b - 4.0 * m * k))*I) / (2.0 * m);
    const double z1 = -k / b;
    bool rst;
    double complex z[1], p[2], zans[1], pans[2];
    c_transfer_function tf;

    // Initialization
    rst = true;
    zans[0] = z1;
    pans[0] = p1;
    pans[1] = p2;
    c_alloc_transfer_function(1, 2, &tf);
    tf.numerator.coefficients[0] = k;
    tf.numerator.coefficients[1] = b;
    tf.denominator.coefficients[0] = k;
    tf.denominator.coefficients[1] = b;
    tf.denominator.coefficients[2] = m;

    // Test
    c_transfer_function_poles(&tf, 2, p);
    c_transfer_function_zeros(&tf, 1, z);

    if (!compare_complex_arrays(1, z, zans, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_poles_zeros -1\n");
    }
    if (!compare_complex_arrays(2, p, pans, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_poles_zeros -2\n");
    }

    // End
    c_free_transfer_function(&tf);
    return rst;
}

bool c_test_tf_multiply()
{
    // Local Variables
    const double tol = 1.0e-8;
    const double m = 1.5;
    const double b = 25.0;
    const double k = 125.3;
    const double r = 2.0e6;
    const double c = 1.5e-6;
    const double x = 3.5;
    bool rst;
    c_transfer_function tf1, tf2, tf;
    double ans1[1], ans2[4];

    // Initialization
    rst = true;
    c_alloc_transfer_function(0, 2, &tf1);
    c_alloc_transfer_function(0, 1, &tf2);
    tf1.numerator.coefficients[0] = 1.0;
    tf1.denominator.coefficients[0] = k;
    tf1.denominator.coefficients[1] = b;
    tf1.denominator.coefficients[2] = m;
    tf2.numerator.coefficients[0] = 1.0;
    tf2.denominator.coefficients[0] = 1.0;
    tf2.denominator.coefficients[1] = r * c;
    ans1[0] = 1.0;
    ans2[0] = k;
    ans2[1] = r * c * k + b;
    ans2[2] = m + r * c * b;
    ans2[3] = r * c * m;

    // Test - Multiplication
    c_transfer_function_multiply(&tf1, &tf2, &tf);
    if (tf.numerator.order != 0)
    {
        rst = false;
        printf("TEST FAILED: c_test_tf_multiply -1\n");
    }
    if (tf.denominator.order != 3)
    {
        rst = false;
        printf("TEST FAILED: c_test_tf_multiply -2\n");
    }
    if (!rst) 
    {
        c_free_transfer_function(&tf1);
        c_free_transfer_function(&tf2);
        c_free_transfer_function(&tf);
        return rst;
    }
    if (!compare_arrays(1, ans1, tf.numerator.coefficients, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_tf_multiply -3\n");
    }
    if (!compare_arrays(4, ans2, tf.denominator.coefficients, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_tf_multiply -4\n");
    }

    // Test - Scaling
    ans1[0] *= x;
    c_scale_transfer_function(x, &tf, &tf);
    if (tf.numerator.order != 0)
    {
        rst = false;
        printf("TEST FAILED: c_test_tf_multiply -5\n");
    }
    if (tf.denominator.order != 3)
    {
        rst = false;
        printf("TEST FAILED: c_test_tf_multiply -6\n");
    }
    if (!rst) 
    {
        c_free_transfer_function(&tf1);
        c_free_transfer_function(&tf2);
        c_free_transfer_function(&tf);
        return rst;
    }
    if (!compare_arrays(1, ans1, tf.numerator.coefficients, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_tf_multiply -7\n");
    }
    if (!compare_arrays(4, ans2, tf.denominator.coefficients, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_tf_multiply -8\n");
    }

    // End
    c_free_transfer_function(&tf1);
    c_free_transfer_function(&tf2);
    c_free_transfer_function(&tf);
    return rst;
}

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

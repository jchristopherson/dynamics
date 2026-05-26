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

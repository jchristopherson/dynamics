#include "dynamics.h"
#include "dynamics_vibrations_tests.h"
#include "dynamics_c_test_helper.h"
#include <math.h>
#include <stdio.h>

bool c_test_q_factor()
{
    // Local Variables
    bool rst;
    double q, ans;
    const double zeta = 0.25;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    ans = 1.0 / (2.0 * zeta);

    // Test
    q = c_q_factor(zeta);
    if (fabs(q - ans) > tol) 
    {
        rst = false;
        printf("TEST FAILED: c_test_q_factor -1\n");
    }

    // End
    return rst;
}


bool c_test_bandwidth()
{
    // Local Variables
    bool rst;
    double df, ans;
    const double fn = 1.0e3;
    const double zeta = 1.0e-1;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    ans = 2.0 * zeta * fn;

    // Test
    df = c_estimate_bandwidth(fn, zeta);
    if (fabs(df - ans) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_bandwidth -1\n");
    }

    // End
    return rst;
}


bool c_test_log_decrement()
{
    // Local Variables
    bool rst;
    double x1, x2, ans, delta;
    const int n = 2;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    x1 = 2.3;
    x2 = 0.8 * x1;
    ans = log(x1 / x2) / n;

    // Test
    delta = c_logarithmic_decrement(x1, x2, n);
    if (fabs(delta - ans) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_log_decrement -1\n");
    }

    // End
    return rst;
}


bool c_test_damping_from_decrement()
{
    // Local Variables;
    bool rst;
    double pi, delta, zeta, ans;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    pi = 2.0 * acos(0.0);
    delta = 0.5;
    ans = delta / sqrt(4.0 * pi * pi + delta * delta);

    // Test
    zeta = c_damping_from_log_decrement(delta);
    if (fabs(zeta - ans) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_damping_from_decrement -1\n");
    }

    // End
    return rst;
}


bool c_test_find_free_rsp_props()
{
    // Local Variables
    bool rst;
    const double fn = 1.0e1;
    const double zeta = 1.0e-1;
    const double xo = 0.0;
    const double vo = 1.0;
    const double dt = 1.0e-3;
    const int n = 1000;
    const double tol = 1.0e-2;
    const double s = 1.0e-2;
    const int np = 2;
    int i;
    double pi, t[n], x[n], wn, wd, A, B, fnx, delta, zx, x1, x2, t1, t2;

    // Initialization
    rst = true;
    pi = 2.0 * acos(0.0);
    wn = 2.0 * pi * fn;
    wd = wn * sqrt(1.0 - zeta * zeta);
    A = (vo + zeta * wn * xo) / wd;
    B = xo;
    for (i = 0; i < n; ++i)
    {
        t[i] = i * dt;
        x[i] = exp(-zeta * wn * t[i]) * (A * sin(wd * t[i]) + B * cos(wd * t[i]));
    }

    // Test
    c_find_free_response_properties(n, t, x, s, np, &delta, &fnx, &x1, &x2, 
        &t1, &t2);
    zx = c_damping_from_log_decrement(delta);

    if (fabs(zx - zeta) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_find_free_rsp_props -1\n");
    }

    if (fabs(fnx - fn) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_find_free_rsp_props -2\n");
    }

    // End
    return rst;
}


bool c_test_rise_time()
{
    // Local Variables
    bool rst;
    const double tol = 1.0e-8;
    const double wn = 1.52e3;
    const double zeta = 0.35;
    double pi, ans, tr, arg;

    // Initialization
    rst = true;
    pi = 2.0 * acos(0.0);
    arg = sqrt(1.0 - zeta * zeta);
    ans = (1.0 / (wn * arg)) * (pi - atan(arg / zeta));

    // Test
    tr = c_rise_time(wn, zeta);
    if (fabs(ans - tr) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_rise_time -1\n");
    }

    // End
    return rst;
}


bool c_test_step_response()
{
    // Local Variables
    bool rst;
    const int n = 1000;
    const double dt = 1.0e-3;
    const double fn = 1.5e1;
    const double zeta = 3.5e-1;
    const double xs = 2.0;
    const double tol = 1.0e-8;
    int i;
    double t[n], x[n], ans[n], wn, wd, A, pi;

    // Initialization
    rst = true;
    wn = 2.0 * pi * fn;
    wd = wn * sqrt(1.0 - zeta * zeta);
    A = zeta * wn / wd;
    for (i = 0; i < n; ++i)
    {
        t[i] = i * dt;
        ans[i] = xs * (1.0 - exp(-zeta * wn * t[i]) * (A * sin(wd * t[i]) + cos(wd * t[i])));
    }

    // Test
    c_evaluate_step_response(n, wn, zeta, xs, t, x);
    if (!compare_arrays(n, ans, x, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_step_response -1\n");
    }

    // End
    return rst;
}


bool c_test_settling_amplitude()
{
    // Local Variables
    bool rst;
    int i;
    const int n = 1000;
    const double tol = 5.0e-2;
    const double dt = 1.0e-3;
    const double fn = 1.3e1;
    const double zeta = 0.25;
    const double xs = 2.0;
    double pi, t[n], x[n], xf, wn;

    // Initialization
    rst = true;
    pi = 2.0 * acos(0.0);
    wn = 2.0 * pi * fn;
    for (i = 0; i < n; ++i)
    {
        t[i] = i * dt;
    }
    c_evaluate_step_response(n, wn, zeta, xs, t, x);

    // Test
    xf = c_find_settling_amplitude(n, x);
    if (fabs(xf - xs) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_settling_amplitude -1\n");
    }

    // End
    return rst;
}


bool c_test_damping_from_overshoot()
{
    // Local Variables
    bool rst;
    const int n = 1000;
    const double tol = 0.1;
    const double dt = 1.0e-3;
    const double fn = 1.2e1;
    const double zeta = 0.82;
    const double xs = 2.0;
    int i;
    double pi, wn, zx, t[n], x[n];

    // Initialization
    rst = true;
    pi = 2.0 * acos(0.0);
    wn = 2.0 * pi * fn;
    for (i = 0; i < n; ++i)
    {
        t[i] = i * dt;
    }
    c_evaluate_step_response(n, wn, zeta, xs, t, x);

    // Test
    zx = c_damping_from_fractional_overshoot(n, x);
    if (fabs(zx - zeta) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_damping_from_overshoot -1\n");
    }

    // End
    return rst;
}

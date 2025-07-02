#include "dynamics.h"
#include "dynamics_frf_tests.h"
#include "dynamics_c_test_helper.h"
#include <math.h>
#include <stdio.h>
#include <float.h>

void c_modal_forcing_term(int n, double freq, double complex *f)
{
    const double complex one = 1.0 + 0.0 * I;
    const double complex zero = 0.0 + 0.0 * I;
    f[0] = one;
    f[1] = zero;
    f[2] = zero;
}

bool c_test_frequency_response()
{
    // Local Variables
    bool rst;
    const int nfreq = 100;
    const double pi = 2.0 * acos(0.0);
    const double tol = 1.0e-3;
    const double fmin = 10.0;
    const double fmax = 1.0e2;
    const double alpha = 0.1;
    const double beta = 2.0e-5;
    const double m1 = 0.5;
    const double m2 = 2.5;
    const double m3 = 0.75;
    const double k1 = 5.0e6;
    const double k2 = 10.0e6;
    const double k3 = 10.0e6;
    const double k4 = 5.0e6;
    int i;
    double df, ans1[3], m[9], k[9], freq[nfreq], omega[nfreq], modes[3], 
        shapes[9];
    double complex rsp[3 * nfreq];
    c_modal_excite fcn;

    // Initialization
    rst = true;
    fcn = c_modal_forcing_term;
    ans1[0] = 2.0 * pi * 232.9225;
    ans1[1] = 2.0 * pi * 749.6189;
    ans1[2] = 2.0 * pi * 923.5669;
    df = (fmax - fmin) / (nfreq - 1.0);
    for (i = 0; i < nfreq; ++i)
    {
        freq[i] = df * i + fmin;
        omega[i] = 2.0 * pi * freq[i];
    }

    // Define the mass matrix
    m[0] = m1;
    m[1] = 0.0;
    m[2] = 0.0;

    m[3] = 0.0;
    m[4] = m2;
    m[5] = 0.0;

    m[6] = 0.0;
    m[7] = 0.0;
    m[8] = m3;

    // Define the stiffness matrix
    k[0] = k1 + k2;
    k[1] = -k2;
    k[2] = 0.0;

    k[3] = -k2;
    k[4] = k2 + k3;
    k[5] = -k3;

    k[6] = 0.0;
    k[7] = -k3;
    k[8] = k3 + k4;

    // Compute the FRF's along with the modal information
    c_frequency_response(3, nfreq, m, 3, k, 3, alpha, beta, omega, fcn, modes, 
        shapes, 3, rsp, nfreq);
    if (!compare_arrays(3, ans1, modes, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_frequency_response -1\n");
    }

    // End
    return rst;
}



bool c_test_modal_response()
{
    // Local Variables
    bool rst;
    double m[9], k[9], ans1[3], vec1[3], vec2[3], vec3[3], freqs[3], shapes[9];
    const double pi = 2.0 * acos(0.0);
    const double m1 = 0.5;
    const double m2 = 2.5;
    const double m3 = 0.75;
    const double k1 = 5.0e6;
    const double k2 = 10.0e6;
    const double k3 = 10.0e6;
    const double k4 = 5.0e6;
    const double tol = 1.0e-3;

    // Initialization
    rst = true;
    ans1[0] = 2.0 * pi * 232.9225;
    ans1[1] = 2.0 * pi * 749.6189;
    ans1[2] = 2.0 * pi * 923.5669;
    
    vec1[0] = 0.7179;
    vec1[1] = 1.0;
    vec1[2] = 0.7466;

    vec2[0] = -0.4192;
    vec2[1] = -0.1638;
    vec2[2] = 1.0;

    vec3[0] = 1.0;
    vec3[1] = -0.1837;
    vec3[2] = 0.1791;

    // Define the mass matrix
    m[0] = m1;
    m[1] = 0.0;
    m[2] = 0.0;

    m[3] = 0.0;
    m[4] = m2;
    m[5] = 0.0;

    m[6] = 0.0;
    m[7] = 0.0;
    m[8] = m3;

    // Define the stiffness matrix
    k[0] = k1 + k2;
    k[1] = -k2;
    k[2] = 0.0;

    k[3] = -k2;
    k[4] = k2 + k3;
    k[5] = -k3;

    k[6] = 0.0;
    k[7] = -k3;
    k[8] = k3 + k4;

    // Compute the modal response
    c_modal_response(3, m, 3, k, 3, freqs, shapes, 3);
    c_normalize_mode_shapes(3, shapes, 3);

    // Tests
    if (!compare_arrays(3, ans1, freqs, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_modal_response -1\n");
    }
    if (!compare_arrays(3, vec1, shapes, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_modal_response -2\n");
    }
    if (!compare_arrays(3, vec2, &shapes[3], tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_modal_response -3\n");
    }
    if (!compare_arrays(3, vec3, &shapes[6], tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_modal_response -4\n");
    }


    // End
    return rst;
}



void c_example_2nd_order_sweep(int n, double freq, double t, const double *x,
    double *dxdt)
{
    // Local Variables
    const double pi = 2.0 * acos(0.0);
    const double z = 1.0e-1;
    const double wn = 2.0 * pi * 50.0;
    double f;

    f = sin(2.0 * pi * freq * t);
    dxdt[0] = x[1];
    dxdt[1] = f - (2.0 * z * wn * x[1] + wn * wn * x[0]);
}

bool c_test_frf_sweep()
{
    // Local Variables
    bool rst;
    const double fmax = 60.0;
    const double fmin = 40.0;
    const double tol = 0.05;
    const int npts = 10;
    const double pi = 2.0 * acos(0.0);
    const double z = 1.0e-1;
    const double wn = 2.0 * pi * 50.0;
    c_frequency_sweep_controls opts;
    int i;
    double df, freq[npts], omega[npts], mag1, mag2, magans1[npts],
        magans2[npts], ref[npts], ratio1[npts], ratio2[npts], init_vals[2];
    double complex rsp[npts * 2], s, tf1, tf2;

    // Initialization
    rst = true;
    c_set_frequency_sweep_defaults(&opts);
    init_vals[0] = 0.0;
    init_vals[1] = 0.0;
    df = (fmax - fmin) / (npts - 1.0);
    for (i = 0; i < npts; ++i)
    {
        freq[i] = df * i + fmin;
        omega[i] = 2.0 * pi * freq[i];
        ref[i] = 1.0;
    }

    // Compute the solution
    for (i = 0; i < npts; ++i)
    {
        s = I * omega[i];
        tf1 = 1.0 / (s * s + 2.0 * z * wn * s + wn * wn);
        tf2 = tf1 * s;
        magans1[i] = cabs(tf1);
        magans2[i] = cabs(tf2);
    }

    // Perform the sweep operation
    c_frf_sweep(2, npts, c_example_2nd_order_sweep, freq, init_vals, 
        DYN_RUNGE_KUTTA_45, rsp, npts, &opts);
    for (i = 0; i < npts; ++i)
    {
        mag1 = cabs(rsp[INDEX(i,0,npts)]);
        mag2 = cabs(rsp[INDEX(i,1,npts)]);
        ratio1[i] = mag1 / magans1[i];
        ratio2[i] = mag2 / magans2[i];
    }

    // Test
    if (!compare_arrays(npts, ref, ratio1, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_frf_sweep -1\n");
    }
    if (!compare_arrays(npts, ref, ratio2, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_frf_sweep -2\n");
    }

    // End
    return rst;
}


bool c_test_frf_fit()
{
    // Local Variables
    bool rst;
    const double pi = 2.0 * acos(0.0);
    const double tol = 5.0e-2;
    const int nfreq = 100;
    const int order = 3;
    const double fmin = 2.0 * pi * 10.0;
    const double fmax = 2.0 * pi * 1.0e3;
    const double alpha = 0.1;
    const double beta = 2.0e-5;
    const double m1 = 0.5;
    const double m2 = 2.5;
    const double m3 = 0.75;
    const double k1 = 5.0e6;
    const double k2 = 10.0e6;
    const double k3 = 10.0e6;
    const double k4 = 5.0e6;
    int i;
    double m[9], k[9], ans[nfreq], mdl[3 * order], fitamp[nfreq], modes[3],
        modeshapes[9], df, freq[nfreq], maxp[3 * order], minp[3 * order],
        maxval;
    double complex fit[nfreq], nrmrsp[3 * nfreq], val;
    c_modal_excite fcn;
    c_iteration_controls controls;
    c_regression_statistics stats[3 * order];

    // Initialization
    rst = true;
    fcn = c_modal_forcing_term;
    df = (fmax - fmin) / (nfreq - 1.0);
    c_set_iteration_controls_defaults(&controls);
    for (i = 0; i < nfreq; ++i)
    {
        freq[i] = fmin + df * i;
    }
    for (i = 0; i < 3 * order; ++i)
    {
        maxp[i] = DBL_MAX;
        minp[i] = -DBL_MAX;
    }

    // Define the mass matrix
    m[0] = m1;
    m[1] = 0.0;
    m[2] = 0.0;

    m[3] = 0.0;
    m[4] = m2;
    m[5] = 0.0;

    m[6] = 0.0;
    m[7] = 0.0;
    m[8] = m3;

    // Define the stiffness matrix
    k[0] = k1 + k2;
    k[1] = -k2;
    k[2] = 0.0;

    k[3] = -k2;
    k[4] = k2 + k3;
    k[5] = -k3;

    k[6] = 0.0;
    k[7] = -k3;
    k[8] = k3 + k4;

    // Compute the frequency response functions
    c_frequency_response(3, nfreq, m, 3, k, 3, alpha, beta, freq, fcn,
        modes, modeshapes, 3, nrmrsp, nfreq);

    // Normalize the response
    val = nrmrsp[0];
    for (i = 0; i < nfreq; ++i)
    {
        nrmrsp[i] /= val;
    }

    // Fit the response
    c_fit_frf(nfreq, order, DYN_RECEPTANCE_MODEL, freq, nrmrsp, maxp, minp,
        &controls, mdl, stats);
    c_evaluate_receptance_frf_model(nfreq, order, mdl, freq, fit);

    // Test
    maxval = 0.0;
    for (i = 0; i < nfreq; ++i)
    {
        fitamp[i] = cabs(fit[i]);
        ans[i] = cabs(nrmrsp[i]);
        maxval = maxval < ans[i] ? ans[i] : maxval;
    }
    if (!compare_arrays(nfreq, fitamp, ans, tol * maxval))
    {
        rst = false;
        printf("TEST FAILED: c_test_frf_fit -1\n");
    }

    // End
    return rst;
}

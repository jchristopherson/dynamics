#include "dynamics.h"
#include "dynamics_c_test_helper.h"
#include "dynamics_kinematics_tests.h"
#include <math.h>
#include <stdio.h>

void inverse_fcn(int njoints, int neqn, const double *x, double *f);

bool c_test_forward_kinematics()
{
    // Local Variables
    bool rst;
    const double L1 = 1.25;
    const double L2 = 5.5;
    const double L3 = 2.0;
    double pi, a1, a2, a3, alpha1, alpha2, alpha3, theta1, theta2, theta3,
        d1, d2, d3, T1[16], T2[16], T3[16], Tref[16], T[16], T23[16];
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    pi = 2.0 * acos(0.0);
    a1 = 0.0;
    a2 = L2;
    a3 = L3;
    alpha1 = 0.5 * pi;
    alpha2 = 0.0;
    alpha3 = 0.0;
    theta1 = 0.0;
    theta2 = 0.5 * pi;
    theta3 = 0.25 * pi;
    d1 = 0.0;
    d2 = -L1;
    d3 = 0.0;

    // Compute the individual transformation matrices
    c_dh_matrix(alpha1, a1, theta1, d1, T1, 4);
    c_dh_matrix(alpha2, a2, theta2, d2, T2, 4);
    c_dh_matrix(alpha3, a3, theta3, d3, T3, 4);

    // Compute the reference matrix
    mtx_mult(4, 4, 4, T2, 4, T3, 4, T23, 4);    // T2 * T3
    mtx_mult(4, 4, 4, T1, 4, T23, 4, Tref, 4);

    // Use the appropriate forward kinematics routine
    c_dh_forward_kinematics_3(T1, 4, T2, 4, T3, 4, T, 4);

    // Test
    if (!compare_matrices(4, 4, Tref, 4, T, 4, tol))
    {
        rst = false;
        printf("TEST FAILED: test_forward_kinematics -1\n");
    }

    // End
    return rst;
}


bool c_test_inverse_kinematics()
{
    // Local Variables
    bool rst;
    const double L1 = 1.25;
    const double L2 = 5.5;
    const double L3 = 2.0;
    double pi, a[3], alpha[3], theta[3], d[3], T[16], ans[3], q[3], qo[3], 
        Tc[16], constraints[6], resid[6];
    const double tol = 1.0e-4;
    c_iteration_behavior ib;
    c_vecfcn fptr;

    // Initialization
    rst = true;
    pi = 2.0 * acos(0.0);
    a[0] = 0.0;
    a[1] = L2;
    a[2] = L3;
    alpha[0] = 0.5 * pi;
    alpha[1] = 0.0;
    alpha[2] = 0.0;
    theta[0] = 0.0;
    theta[1] = 0.5 * pi;
    theta[2] = 0.25 * pi;
    d[0] = 0.0;
    d[1] = -L1;
    d[2] = 0.0;
    ans[0] = theta[0];
    ans[1] = theta[1];
    ans[2] = theta[2];

    // Solve the forward kinematics problem to define the constraints
    c_dh_forward_kinematics(3, alpha, a, theta, d, T, 4);
    constraints[0] = T[INDEX(0,3,4)];
    constraints[1] = T[INDEX(1,3,4)];
    constraints[2] = T[INDEX(2,3,4)];
    constraints[3] = T[INDEX(0,0,4)];
    constraints[4] = T[INDEX(1,0,4)];
    constraints[5] = T[INDEX(2,0,4)];

    // Define an initial guess at the joint variables
    qo[0] = 0.0;
    qo[1] = 0.0;
    qo[2] = 0.0;

    // Solver the inverse problem
    fptr = inverse_fcn;
    c_solve_inverse_kinematics(3, 6, fptr, qo, constraints, q, resid, &ib);

    // Test
    if (!compare_arrays(3, ans, q, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_inverse_kinematics -1\n");
    }

    // End
    return rst;
}

void inverse_fcn(int njoints, int neqn, const double *x, double *f)
{
    // Local Variables
    const double L1 = 1.25;
    const double L2 = 5.5;
    const double L3 = 2.0;
    double pi, T[16], alpha[3], a[3], d[3];

    // Initialization
    pi = 2.0 * acos(0.0);
    a[0] = 0.0;
    a[1] = L2;
    a[2] = L3;
    alpha[0] = 0.5 * pi;
    alpha[1] = 0.0;
    alpha[2] = 0.0;
    d[0] = 0.0;
    d[1] = -L1;
    d[2] = 0.0;

    // Define the location of the end-effector
    c_dh_forward_kinematics(3, alpha, a, x, d, T, 4);

    // Define the position of the end-effector
    f[0] = T[INDEX(0,3,4)];
    f[1] = T[INDEX(1,3,4)];
    f[2] = T[INDEX(2,3,4)];

    // Define it's orientation
    f[3] = T[INDEX(0,0,4)];
    f[4] = T[INDEX(1,0,4)];
    f[5] = T[INDEX(2,0,4)];
}


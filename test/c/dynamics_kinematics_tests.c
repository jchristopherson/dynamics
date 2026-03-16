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

bool c_test_define_link_csys()
{
    bool rst;
    double xim1[3], zim1[3], zi[3], rim1[3], ri[3], u[3], ians[3], jans[3],
        kans[3];
    c_coordinate_system csys;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    xim1[0] = 1.0;  xim1[1] = 0.0;  xim1[2] = 0.0;
    zim1[0] = 0.0;  zim1[1] = 0.0;  zim1[2] = 1.0;
    zi[0] = 1.0;    zi[1] = 0.0;    zi[2] = 0.0;
    rim1[0] = 0.0;  rim1[1] = 0.0;  rim1[2] = 0.0;
    ri[0] = 1.0;    ri[1] = 1.0;    ri[2] = 0.0;
    u[0] = 0.0;     u[1] = 1.0;     u[2] = 0.0;
    ians[0] = 0.0;  ians[1] = 1.0;  ians[2] = 0.0;
    kans[0] = 1.0;  kans[1] = 0.0;  kans[2] = 0.0;
    c_cross_product(kans, ians, jans);

    // Build the coordinate system
    c_define_link_csys(xim1, zim1, zi, rim1, ri, &csys);

    // Tests
    if (!compare_arrays(3, ians, csys.i, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_define_link_csys -1\n");
    }

    if (!compare_arrays(3, jans, csys.j, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_define_link_csys -2\n");
    }

    if (!compare_arrays(3, kans, csys.k, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_define_link_csys -3\n");
    }

    if (!compare_arrays(3, u, csys.origin, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_define_link_csys -4\n");
    }

    // End
    return rst;
}

bool c_test_define_csys()
{
    bool rst;
    double i[3], j[3], k[3], o[3];
    c_coordinate_system csys;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    i[0] = 1.0;     i[1] = 0.0;     i[2] = 0.0;
    j[0] = 0.0;     j[1] = 1.0;     j[2] = 0.0;
    k[0] = 0.0;     k[1] = 0.0;     k[2] = 0.0;
    create_random_vector(3, o);
    c_define_csys(i, j, k, o, &csys);

    // Tests
    if (!compare_arrays(3, i, csys.i, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_define_csys -1\n");
    }

    if (!compare_arrays(3, j, csys.j, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_define_csys -2\n");
    }

    if (!compare_arrays(3, k, csys.k, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_define_csys -3\n");
    }

    if (!compare_arrays(3, o, csys.origin, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_define_csys -4\n");
    }

    // End
    return rst;
}

bool c_test_build_dh_table()
{
    bool rst;
    const double tol = 1.0e-8;

    // Local Variables
    const double L1 = 1.75;
    const double L2 = 3.5;
    const double L3 = 2.0;
    const double theta1 = 0.125;
    const double theta2 = -0.25;
    const double theta3 = 1.25;
    c_coordinate_system csys[4];
    c_dh_table tbl;
    int ii;
    double o[3], ob1[3], ob2[3], ob3[3], o1[3], o2[3], o3[3], ib[3], jb[3], 
        kb[3], i1[3], j1[3], i2[3], j2[3], i3[3], j3[3], R1[9], R2[9], R3[9], 
        R12[9], R123[9], T[16], p3b[4], p3[4], pans[4];

    // Initialization
    rst = true;

    // Define the origin
    o[0] = 0.0;     o[1] = 0.0;     o[2] = 0.0;

    // Define the matrices used to locate the joints
    c_rotate_z(theta1, R1, 3);
    c_rotate_z(theta2, R2, 3);
    c_rotate_z(theta3, R3, 3);
    c_matmul(3, 3, 3, 1.0, R1, 3, R2, 3, 0.0, R12, 3);  // = R1 * R2
    c_matmul(3, 3, 3, 1.0, R12, 3, R3, 3, 0.0, R123, 3); // = R1 * R2 * R3

    // Define the locations of each of the joints
    ob1[0] = L1;    ob1[1] = 0.0;   ob1[2] = 0.0;
    ob2[0] = L2;    ob2[1] = 0.0;   ob2[2] = 0.0;
    ob3[0] = L3;    ob3[1] = 0.0;   ob3[2] = 0.0;
    p3b[0] = 0.0;   p3b[1] = 0.0;   p3b[2] = 0.0;   p3b[3] = 1.0;

    c_matmul(3, 1, 3, 1.0, R1, 3, ob1, 3, 0.0, o1, 3);   // o1 = R1 * ob1
    c_matmul(3, 1, 3, 1.0, R12, 3, ob2, 3, 0.0, o2, 3);  // o2 = R1 * R2 * ob2
    c_matmul(3, 1, 3, 1.0, R123, 3, ob3, 3, 0.0, o3, 3); // o3 = R1 * R2 * R3 * ob3

    for (ii = 0; ii < 3; ++ii)
    {
        o2[ii] += o1[ii];
        o3[ii] += o2[ii];
    }
    pans[0] = o3[0];    pans[1] = o3[1];    pans[2] = o3[2];    pans[3] = 1.0;

    // Define unit vectors
    ib[0] = 1.0;    ib[1] = 0.0;    ib[2] = 0.0;
    jb[0] = 0.0;    jb[1] = 1.0;    jb[2] = 0.0;
    kb[0] = 0.0;    kb[1] = 0.0;    kb[2] = 1.0;

    // Define the ground coordinate system
    c_define_csys(ib, jb, kb, o, &csys[0]);

    // Define CSYS 1
    c_matmul(3, 1, 3, 1.0, R1, 3, ib, 3, 0.0, i1, 3);
    c_matmul(3, 1, 3, 1.0, R1, 3, jb, 3, 0.0, j1, 3);
    c_define_csys(i1, j1, kb, o1, &csys[1]);
    
    // Define CSYS 2
    c_matmul(3, 1, 3, 1.0, R12, 3, ib, 3, 0.0, i2, 3);
    c_matmul(3, 1, 3, 1.0, R12, 3, jb, 3, 0.0, j2, 3);
    c_define_csys(i2, j2, kb, o2, &csys[2]);

    // Define CSYS 3
    c_matmul(3, 1, 3, 1.0, R123, 3, ib, 3, 0.0, i3, 3);
    c_matmul(3, 1, 3, 1.0, R123, 3, jb, 3, 0.0, j3, 3);
    c_define_csys(i3, j3, kb, o3, &csys[3]);

    // Build the table
    c_build_dh_table(4, csys, &tbl);

    // Build the linkage transformation matrix & locate the end effector
    c_dh_forward_kinematics_table(&tbl, T, 4);
    c_matmul(4, 1, 4, 1.0, T, 4, p3b, 4, 0.0, p3, 4);

    // Tests
    if (tbl.count != 3)
    {
        printf("TEST FAILED: c_test_build_dh_table\n");
        c_free_dh_table(&tbl);
        return false;
    }
    
    if (fabs(tbl.parameters[0].link_length - L1) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_build_dh_table -1a\n");
    }
    if (fabs(tbl.parameters[0].link_twist) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_build_dh_table -1b\n");
    }
    if (fabs(tbl.parameters[0].link_offset) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_build_dh_table -1c\n");
    }
    if (fabs(tbl.parameters[0].joint_angle - theta1) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_build_dh_table -1d\n");
    }

    if (fabs(tbl.parameters[1].link_length - L2) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_build_dh_table -2a\n");
    }
    if (fabs(tbl.parameters[1].link_twist) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_build_dh_table -2b\n");
    }
    if (fabs(tbl.parameters[1].link_offset) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_build_dh_table -2c\n");
    }
    if (fabs(tbl.parameters[1].joint_angle - theta2) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_build_dh_table -2d\n");
    }

    if (fabs(tbl.parameters[2].link_length - L3) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_build_dh_table -3a\n");
    }
    if (fabs(tbl.parameters[2].link_twist) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_build_dh_table -3b\n");
    }
    if (fabs(tbl.parameters[2].link_offset) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_build_dh_table -3c\n");
    }
    if (fabs(tbl.parameters[2].joint_angle - theta3) > tol)
    {
        rst = false;
        printf("TEST FAILED: c_test_build_dh_table -3d\n");
    }

    if (!compare_arrays(4, p3, pans, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_build_dh_table -4\n");
    }

    // Clean up
    c_free_dh_table(&tbl);

    // End
    return rst;
}


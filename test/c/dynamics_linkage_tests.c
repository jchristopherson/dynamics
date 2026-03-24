#include "dynamics_c_test_helper.h"
#include "dynamics_linkage_tests.h"
#include <stdio.h>
#include <math.h>

c_serial_linkage build_serial_linkage_1()
{
    const double pi = 2.0 * acos(0.0);
    const double Lo = 1.5;
    const double dstatic = 1.25;
    c_serial_linkage linkage;
    c_binary_link links[3];
    int i, j, k;

    // Link 1
    links[0].link_twist = -0.5 * pi;
    links[0].link_offset = Lo;
    links[0].joint_angle = 0.0;
    links[0].link_length = 0.0;
    links[0].joint_type = DYN_REVOLUTE_JOINT;
    links[0].mass = 1.0;

    // Link 2
    links[1].link_twist = 0.5 * pi;
    links[1].link_offset = 0.0;
    links[1].joint_angle = 0.0;
    links[1].link_length = 0.0;
    links[1].joint_type = DYN_REVOLUTE_JOINT;
    links[1].mass = 1.0;

    // Link 3
    links[2].link_twist = 0.0;
    links[2].link_offset = dstatic;
    links[2].joint_angle = 0.0;
    links[2].link_length = 0.0;
    links[2].joint_type = DYN_PRISMATIC_JOINT;
    links[2].mass = 1.0;

    // Define the CG and inertia terms
    for (k = 0; k < 3; ++k)
    {
        for (i = 0; i < 3; ++i)
        {
            links[k].cg[i] = 0.0;
        }
        for (j = 0; j < 3; ++j)
        {
            for (i = 0; i < 3; ++i)
            {
                links[k].inertia[3 * j + i] = i == j ? 1.0 : 0.0;
            }
        }
    }

    // Construct the linkage
    c_build_serial_linkage(3, links, &linkage);

    // End
    return linkage;
}

bool c_test_serial_linkage_forward_kinematics()
{
    bool rst;
    int i, n;
    double q[3], d[3], alpha[3], a[3], theta[3], T1[16], T2[16], T3[16], 
        Tans[16], T[16];
    c_serial_linkage linkage;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    linkage = build_serial_linkage_1();
    n = linkage.link_count;
    create_random_vector(3, q);
    for (i = 0; i < n; ++i)
    {
        alpha[i] = linkage.links[i].link_twist;
        a[i] = linkage.links[i].link_length;
        if (linkage.links[i].joint_type == DYN_PRISMATIC_JOINT)
        {
            theta[i] = linkage.links[i].joint_angle;
            d[i] = linkage.links[i].link_offset + q[i];
        }
        else
        {
            theta[i] = linkage.links[i].joint_angle + q[i];
            d[i] = linkage.links[i].link_offset;
        }
    }

    // Compute the solution
    c_dh_matrix(alpha[0], a[0], theta[0], d[0], T1, 4);
    c_dh_matrix(alpha[1], a[1], theta[1], d[1], T2, 4);
    c_dh_matrix(alpha[2], a[2], theta[2], d[2], T3, 4);
    c_dh_forward_kinematics_3(T1, 4, T2, 4, T3, 4, Tans, 4);

    // Evaluate the linkage forward kinematics
    c_serial_linkage_forward_kinematics(3, &linkage, q, T, 4);

    // Test
    if (!compare_matrices(4, 4, Tans, 4, T, 4, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_serial_linkage_forward_kinematics -1\n");
    }

    // End
    c_free_serial_linkage(&linkage);
    return rst;
}

bool c_test_serial_linkage_jacobian()
{
    bool rst;
    int i, n, jtypes[3];
    double q[3], d[3], alpha[3], a[3], theta[3], J[18], Jans[18];
    c_serial_linkage linkage;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    linkage = build_serial_linkage_1();
    n = linkage.link_count;
    create_random_vector(3, q);
    for (i = 0; i < n; ++i)
    {
        alpha[i] = linkage.links[i].link_twist;
        a[i] = linkage.links[i].link_length;
        if (linkage.links[i].joint_type == DYN_PRISMATIC_JOINT)
        {
            theta[i] = linkage.links[i].joint_angle;
            d[i] = linkage.links[i].link_offset + q[i];
            jtypes[i] = DYN_PRISMATIC_JOINT;
        }
        else
        {
            theta[i] = linkage.links[i].joint_angle + q[i];
            d[i] = linkage.links[i].link_offset;
            jtypes[i] = DYN_REVOLUTE_JOINT;
        }
    }

    // Compute the solution
    c_dh_jacobian(n, alpha, a, theta, d, jtypes, Jans, 6);
    
    // Evaluate the linkage jacobian
    c_serial_linkage_jacobian(n, &linkage, q, J, 6);

    // Test
    if (!compare_matrices(6, n, Jans, 6, J, 6, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_serial_linkage_jacobian -1\n");
    }

    // End
    c_free_serial_linkage(&linkage);
    return rst;
}

bool c_test_serial_linkage_inverse_kinematics()
{
    bool rst;
    int i, n;
    double qo[3], q[3], s[3], d[3], alpha[3], a[3], theta[3], T1[16], T2[16], 
        T3[16], T[16];
    c_serial_linkage linkage;
    c_iteration_behavior ib;
    const double tol = 1.0e-8;

    // Initialization
    rst = true;
    linkage = build_serial_linkage_1();
    n = linkage.link_count;
    create_random_vector(3, q); // also the solution vector
    for (i = 0; i < n; ++i)
    {
        qo[i] = 0.0;
        alpha[i] = linkage.links[i].link_twist;
        a[i] = linkage.links[i].link_length;
        if (linkage.links[i].joint_type == DYN_PRISMATIC_JOINT)
        {
            theta[i] = linkage.links[i].joint_angle;
            d[i] = linkage.links[i].link_offset + q[i];
        }
        else
        {
            theta[i] = linkage.links[i].joint_angle + q[i];
            d[i] = linkage.links[i].link_offset;
        }
    }

    // Compute the target
    c_dh_matrix(alpha[0], a[0], theta[0], d[0], T1, 4);
    c_dh_matrix(alpha[1], a[1], theta[1], d[1], T2, 4);
    c_dh_matrix(alpha[2], a[2], theta[2], d[2], T3, 4);
    c_dh_forward_kinematics_3(T1, 4, T2, 4, T3, 4, T, 4);

    // Solve the inverse model
    c_serial_linkage_inverse_kinematics(n, &linkage, qo, T, 4, s, &ib);

    // Test
    if (!compare_arrays(n, q, s, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_serial_linkage_inverse_kinematics -1\n");
    }

    // End
    c_free_serial_linkage(&linkage);
    return rst;
}

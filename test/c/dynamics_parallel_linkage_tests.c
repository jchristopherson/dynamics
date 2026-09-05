#include "dynamics_c_test_helper.h"
#include "dynamics_parallel_linkage_tests.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define PI 3.14159265358979323846

static void identity(double *T)
{
    int i;
    for (i = 0; i < 16; ++i) T[i] = 0.0;
    T[0] = 1.0;
    T[5] = 1.0;
    T[10] = 1.0;
    T[15] = 1.0;
}

static int make_link(double length, c_mechanism_link *link)
{
    if (c_alloc_mechanism_link(2, link) != 0) return -1;
    identity(&link->frames[0]);
    identity(&link->frames[16]);
    link->frames[28] = length;
    return 0;
}

static void make_four_bar(c_mechanism_link links[4], c_joint joints[4])
{
    make_link(4.0, &links[0]);
    make_link(1.0, &links[1]);
    make_link(3.5, &links[2]);
    make_link(3.0, &links[3]);

    joints[0].joint_type = DYN_REVOLUTE_JOINT;
    joints[0].parent_link = 1;
    joints[0].parent_frame = 1;
    joints[0].child_link = 2;
    joints[0].child_frame = 1;
    joints[0].actuated = true;

    joints[1].joint_type = DYN_REVOLUTE_JOINT;
    joints[1].parent_link = 2;
    joints[1].parent_frame = 2;
    joints[1].child_link = 3;
    joints[1].child_frame = 1;
    joints[1].actuated = false;

    joints[2].joint_type = DYN_REVOLUTE_JOINT;
    joints[2].parent_link = 3;
    joints[2].parent_frame = 2;
    joints[2].child_link = 4;
    joints[2].child_frame = 2;
    joints[2].actuated = false;

    joints[3].joint_type = DYN_REVOLUTE_JOINT;
    joints[3].parent_link = 4;
    joints[3].parent_frame = 1;
    joints[3].child_link = 1;
    joints[3].child_frame = 2;
    joints[3].actuated = false;
}

static void four_bar_configuration(double theta, double q[4], double point[2])
{
    const double coupler = 3.5;
    const double rocker = 3.0;
    const double ground = 4.0;
    double b[2], d[2], u[2], w[2], length, x, height, psi, phi;

    b[0] = cos(theta);
    b[1] = sin(theta);
    d[0] = ground;
    d[1] = 0.0;
    length = hypot(d[0] - b[0], d[1] - b[1]);
    u[0] = (d[0] - b[0]) / length;
    u[1] = (d[1] - b[1]) / length;
    w[0] = -u[1];
    w[1] = u[0];
    x = 0.5 * (length * length + coupler * coupler - rocker * rocker) /
        length;
    height = sqrt(coupler * coupler - x * x);
    point[0] = b[0] + x * u[0] + height * w[0];
    point[1] = b[1] + x * u[1] + height * w[1];
    psi = atan2(point[1] - b[1], point[0] - b[0]);
    phi = atan2(point[1], point[0] - ground);
    q[0] = theta;
    q[1] = psi - theta;
    q[2] = phi - psi;
    q[3] = -phi;
}

bool c_test_parallel_linkage_topology()
{
    bool rst = true;
    c_mechanism_link links[4];
    c_joint joints[4];
    c_mechanism linkage;
    double tool[16];

    make_four_bar(links, joints);
    identity(tool);
    linkage = c_create_parallel_linkage(4, links, 4, joints, 1, 3, tool, 4);
    if (c_mechanism_link_count(linkage) != 4 ||
        c_mechanism_joint_count(linkage) != 4 ||
        c_mechanism_variable_count(linkage) != 4 ||
        c_mechanism_loop_count(linkage) != 1 ||
        c_mechanism_constraint_count(linkage) != 6 ||
        c_mechanism_space_dimension(linkage) != 6)
    {
        rst = false;
        printf("TEST FAILED: c_test_parallel_linkage_topology -1\n");
    }

    c_free_mechanism(linkage);
    for (int i = 0; i < 4; ++i) c_free_mechanism_link(&links[i]);
    return rst;
}

bool c_test_planar_linkage_kinematics()
{
    bool rst = true;
    const double tol = 1.0e-8;
    const double theta = 0.6;
    c_mechanism_link links[4];
    c_joint joints[4];
    c_mechanism linkage;
    c_iteration_behavior ib;
    double tool[16], q[4], qactual[4], point[2], T[16], f[3];
    double constraint_jacobian[12], jacobian[3];

    make_four_bar(links, joints);
    identity(tool);
    tool[12] = 3.5;
    linkage = c_create_planar_linkage(4, links, 4, joints, 1, 3, tool, 4);

    if (c_mechanism_space_dimension(linkage) != 3 ||
        c_mechanism_degrees_of_freedom(linkage) != 1 ||
        c_mechanism_actuated_variable_count(linkage) != 1 ||
        c_mechanism_link_frame_count(linkage, 2) != 2)
    {
        rst = false;
        printf("TEST FAILED: c_test_planar_linkage_kinematics -1\n");
    }

    four_bar_configuration(theta, q, point);
    c_mechanism_set_configuration(linkage, 4, q);
    c_mechanism_get_configuration(linkage, 4, qactual);
    c_mechanism_constraints(linkage, 4, q, 3, f);
    if (!compare_arrays(4, q, qactual, tol) ||
        !compare_arrays(3, f, (double[3]){0.0, 0.0, 0.0}, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_planar_linkage_kinematics -2\n");
    }

    c_mechanism_forward_kinematics(linkage, 1, &theta, T, 4, &ib);
    if (fabs(T[12] - point[0]) > 1.0e-6 || fabs(T[13] - point[1]) > 1.0e-6)
    {
        rst = false;
        printf("TEST FAILED: c_test_planar_linkage_kinematics -3\n");
    }

    c_mechanism_constraint_jacobian(linkage, 4, q, constraint_jacobian, 3);
    c_mechanism_jacobian(linkage, 1, &theta, jacobian, 3);
    c_mechanism_solve_configuration(linkage, 1, &theta, 4, qactual, &ib);
    if (!compare_arrays(1, qactual, &theta, 1.0e-6))
    {
        rst = false;
        printf("TEST FAILED: c_test_planar_linkage_kinematics -4\n");
    }

    c_free_mechanism(linkage);
    for (int i = 0; i < 4; ++i) c_free_mechanism_link(&links[i]);
    return rst;
}
#include "dynamics_variational_tests.h"
#include "dynamics.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static void fixed_force(const c_variational_state *state, double *force,
    double *torque, void *user_data)
{
    (void)state; (void)user_data;
    force[0] = 0.0; force[1] = 0.0; force[2] = -19.62;
    torque[0] = 0.0; torque[1] = 0.0; torque[2] = 0.0;
}

static void fixed_constraint(const c_variational_state *state, int n,
    double *value, void *user_data)
{
    const double *position = state->position;
    (void)user_data;
    for (int i = 0; i < n; ++i) value[i] = position[i];
}

static void fixed_jacobian(const c_variational_state *state, int n,
    double *jacobian, int ldj, void *user_data)
{
    (void)state; (void)user_data;
    for (int i = 0; i < ldj * 6; ++i) jacobian[i] = 0.0;
    for (int i = 0; i < n; ++i) jacobian[i * ldj + i] = 1.0;
}

bool c_test_variational_integrator(void)
{
    c_rigid_body body = {0};
    c_variational_integrator_settings settings;
    c_quaternion q0 = {1.0, 0.0, 0.0, 0.0};
    c_quaternion q[3];
    double p0[3] = {0.0, 0.0, 0.0};
    double v0[3] = {0.0, 0.0, 0.0};
    double w0[3] = {0.0, 0.0, 0.0};
    double p[9], v[9], w[9], multipliers[9];

    body.mass = 2.0;
    body.inertia[0] = body.inertia[4] = body.inertia[8] = 1.0;
    c_default_variational_integrator_settings(&settings);
    c_variational_integrator_solve(1, &body, 3, 0.01, p0, &q0, v0, w0, 3,
        fixed_force, fixed_constraint, fixed_jacobian, NULL, &settings,
        p, q, v, w, multipliers);
    if (fabs(p[6]) > 1.0e-10 || fabs(p[7]) > 1.0e-10 ||
        fabs(p[8]) > 1.0e-10 || fabs(multipliers[2] - 19.62) > 1.0e-6)
    {
        printf("TEST FAILED: c_test_variational_integrator\n");
        return false;
    }
    return true;
}

static void identity(double *matrix)
{
    for (int i = 0; i < 16; ++i) matrix[i] = 0.0;
    matrix[0] = matrix[5] = matrix[10] = matrix[15] = 1.0;
}

static int make_massive_link(double length, double mass, c_mechanism_link *link)
{
    if (c_alloc_mechanism_link(2, link) != 0) return -1;
    identity(link->frames);
    identity(link->frames + 16);
    link->frames[28] = length;
    link->mass = mass;
    link->cg[0] = 0.5 * length;
    link->inertia[0] = 1.0e-3;
    link->inertia[4] = link->inertia[8] = mass * length * length / 12.0;
    return 0;
}

static void four_bar_configuration(double theta, double q[4])
{
    const double coupler = 3.5, rocker = 3.0, ground = 4.0;
    double b[2] = {cos(theta), sin(theta)};
    double d[2] = {ground, 0.0};
    double length = hypot(d[0] - b[0], d[1] - b[1]);
    double u[2] = {(d[0] - b[0]) / length, (d[1] - b[1]) / length};
    double x = 0.5 * (length * length + coupler * coupler - rocker * rocker) /
        length;
    double height = sqrt(coupler * coupler - x * x);
    double point[2] = {b[0] + x*u[0] - height*u[1],
        b[1] + x*u[1] + height*u[0]};
    double psi = atan2(point[1] - b[1], point[0] - b[0]);
    double phi = atan2(point[1], point[0] - ground);
    q[0] = theta; q[1] = psi - theta; q[2] = phi - psi; q[3] = -phi;
}

static double fixed_motion(double t, void *user_data)
{
    (void)t;
    return *(const double*)user_data;
}

bool c_test_linkage_dynamics(void)
{
    c_mechanism_link links[4];
    c_joint joints[4] = {0};
    c_variational_integrator_settings settings;
    c_linkage_dynamic_model model;
    c_linkage_dynamic_model serial_model;
    c_serial_linkage serial = {0};
    c_joint_reaction reactions[4];
    double configuration[4], theta = 0.7;
    double gravity[3] = {0.0, -9.81, 0.0};
    double force[9] = {0.0}, torque[9] = {0.0};
    double position[18], velocity[18], angular_velocity[18];
    c_quaternion orientation[6];
    double multipliers[36]; /* (17 linkage + 1 prescribed) by 2 times */
    bool result = true;

    c_alloc_serial_linkage(1, &serial);
    serial.links[0].joint_type = DYN_REVOLUTE_JOINT;
    serial.links[0].link_length = 1.0;
    serial.links[0].mass = 1.0;
    serial.links[0].cg[0] = -0.5;
    serial.links[0].inertia[0] = serial.links[0].inertia[4] =
        serial.links[0].inertia[8] = 1.0;
    serial_model = c_create_serial_linkage_dynamic_model(&serial, 1, &theta);
    if (!serial_model || c_linkage_dynamic_body_count(serial_model) != 1 ||
        c_linkage_dynamic_joint_count(serial_model) != 1)
        result = false;
    c_free_linkage_dynamic_model(serial_model);
    c_free_serial_linkage(&serial);

    make_massive_link(4.0, 2.0, &links[0]);
    make_massive_link(1.0, 1.0, &links[1]);
    make_massive_link(3.5, 1.5, &links[2]);
    make_massive_link(3.0, 1.2, &links[3]);
    int pairs[4][2] = {{1,2},{2,3},{3,4},{4,1}};
    int frames[4][2] = {{1,1},{2,1},{2,2},{1,2}};
    for (int i = 0; i < 4; ++i) {
        joints[i].joint_type = DYN_REVOLUTE_JOINT;
        joints[i].parent_link = pairs[i][0]; joints[i].child_link = pairs[i][1];
        joints[i].parent_frame = frames[i][0]; joints[i].child_frame = frames[i][1];
    }
    four_bar_configuration(theta, configuration);
    model = c_create_linkage_dynamic_model(true, 4, links, 4, joints, 1, 4,
        configuration);
    c_default_variational_integrator_settings(&settings);
    if (!model || c_linkage_dynamic_body_count(model) != 3 ||
        c_linkage_dynamic_joint_count(model) != 4 ||
        c_linkage_dynamic_constraint_count(model) != 17) result = false;
    if (result) {
        c_linkage_dynamic_solve(model, 3, 18, &settings, 2, 1.0e-4, gravity,
            force, torque, 1, fixed_motion, &theta, position, orientation,
            velocity, angular_velocity, multipliers);
        c_linkage_dynamic_joint_reactions(model, 3, 18, 1.0e-4,
            position + 9, orientation + 3, velocity + 9,
            angular_velocity + 9, multipliers + 18, reactions);
        for (int i = 0; i < 3; ++i)
            if (!isfinite(reactions[0].force[i])) result = false;
    }
    if (!result) printf("TEST FAILED: c_test_linkage_dynamics\n");
    c_free_linkage_dynamic_model(model);
    for (int i = 0; i < 4; ++i) c_free_mechanism_link(&links[i]);
    return result;
}

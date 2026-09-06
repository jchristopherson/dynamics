#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <dynamics.h>

/*

This example analyzes a planar four-bar linkage.  The mechanism contains a
single closed kinematic loop, so its forward kinematics require the solution of
the loop-closure constraints rather than a simple accumulation of link
transformations.

This is the C equivalent of the four_bar_example_1 Fortran example, though the
results are written to the terminal rather than plotted.

*/

#define NPTS 13
#define PI 3.14159265358979323846

// Populates a 4-by-4 identity matrix stored in column-major format.
static void identity(double *T)
{
    int i;
    for (i = 0; i < 16; ++i) T[i] = 0.0;
    T[0] = 1.0;
    T[5] = 1.0;
    T[10] = 1.0;
    T[15] = 1.0;
}

// Defines a planar link of the requested length carrying a joint frame at each
// end.  The body frame is coincident with the first of the two frames.
static int planar_link(double length, c_mechanism_link *lnk)
{
    if (c_alloc_mechanism_link(2, lnk) != 0) return -1;
    identity(&lnk->frames[0]);
    identity(&lnk->frames[16]);
    lnk->frames[16 + 12] = length;  // the x-coordinate of the second frame
    return 0;
}

int main()
{
    // Model Properties
    const double crank = 1.0;
    const double coupler = 3.5;
    const double rocker = 3.0;
    const double ground = 4.0;

    // Local Variables
    int i, j, k, nlinks, nvar, nframes;
    double theta, tool[16], T[16], F[16], P[16], jac[3];
    double *q;
    c_mechanism_link links[4];
    c_joint joints[4];
    c_mechanism linkage;
    c_iteration_behavior ib;
    double config[4] = {0.0, 0.5 * PI, -0.5 * PI, 0.0};

    // Define the links
    planar_link(ground, &links[0]);
    planar_link(crank, &links[1]);
    planar_link(coupler, &links[2]);
    planar_link(rocker, &links[3]);

    // Connect the links.  The crank is the driven member.
    joints[0].joint_type = DYN_REVOLUTE_JOINT;
    joints[0].parent_link = 1;  joints[0].parent_frame = 1;
    joints[0].child_link = 2;   joints[0].child_frame = 1;
    joints[0].actuated = true;

    joints[1].joint_type = DYN_REVOLUTE_JOINT;
    joints[1].parent_link = 2;  joints[1].parent_frame = 2;
    joints[1].child_link = 3;   joints[1].child_frame = 1;
    joints[1].actuated = false;

    joints[2].joint_type = DYN_REVOLUTE_JOINT;
    joints[2].parent_link = 3;  joints[2].parent_frame = 2;
    joints[2].child_link = 4;   joints[2].child_frame = 2;
    joints[2].actuated = false;

    joints[3].joint_type = DYN_REVOLUTE_JOINT;
    joints[3].parent_link = 4;  joints[3].parent_frame = 1;
    joints[3].child_link = 1;   joints[3].child_frame = 2;
    joints[3].actuated = false;

    // The end-effector is the distal end of the coupler
    identity(tool);
    tool[12] = coupler;

    // Build the mechanism
    linkage = c_create_planar_linkage(4, links, 4, joints, 1, 3, tool, 4);

    // Display the mobility of the mechanism
    nlinks = c_mechanism_link_count(linkage);
    nvar = c_mechanism_variable_count(linkage);
    printf("Number of independent loops: %i\n", c_mechanism_loop_count(linkage));
    printf("Number of joint variables: %i\n", nvar);
    printf("Number of constraint equations: %i\n",
        c_mechanism_constraint_count(linkage));
    printf("Degrees of freedom: %i\n",
        c_mechanism_degrees_of_freedom(linkage));

    // Establish a starting configuration.  A closed-loop mechanism admits more
    // than one assembly mode, so the starting estimate selects the branch.
    c_mechanism_set_configuration(linkage, 4, config);

    // Sweep the crank and trace the coupler point
    printf("\n  CRANK          X          Y      ANGLE\n");
    for (i = 0; i < NPTS; ++i)
    {
        theta = 2.0 * PI * (double)i / (double)(NPTS - 1);
        c_mechanism_forward_kinematics(linkage, 1, &theta, T, 4, &ib);
        printf("%7.3f %10.4f %10.4f %10.4f\n", theta, T[12], T[13],
            atan2(T[1], T[0]));
    }

    // Examine the sensitivity of the coupler point to the crank at mid-stroke
    theta = 0.5 * PI;
    c_mechanism_jacobian(linkage, 1, &theta, jac, 3);
    printf("\nJacobian at a crank angle of 90 degrees:\n");
    for (i = 0; i < 3; ++i) printf("%20.15f\n", jac[i]);

    // Locate each joint with the crank at 45 degrees.  The complete set of
    // joint variables, not just the actuated variable, is required in order to
    // locate each link.
    theta = 0.25 * PI;
    q = (double*)malloc((size_t)nvar * sizeof(double));
    c_mechanism_solve_configuration(linkage, 1, &theta, nvar, q, &ib);

    printf("\nLinkage geometry with the crank at 45 degrees:\n");
    for (j = 1; j <= nlinks; ++j)
    {
        // Locate the link, then locate each joint frame the link carries
        c_mechanism_body_transform(linkage, j, nvar, q, T, 4);
        nframes = c_mechanism_link_frame_count(linkage, j);
        printf("Link %i:", j);
        for (k = 1; k <= nframes; ++k)
        {
            c_mechanism_link_frame(linkage, j, k, F, 4);
            c_matmul(4, 4, 4, 1.0, T, 4, F, 4, 0.0, P, 4);
            printf("  (%8.4f, %8.4f)", P[12], P[13]);
        }
        printf("\n");
    }

    // Verify that the loop-closure constraints are satisfied
    c_mechanism_constraints(linkage, nvar, q, 3, jac);
    printf("\nLoop-closure residual: %.3e\n",
        fmax(fmax(fabs(jac[0]), fabs(jac[1])), fabs(jac[2])));

    // Clean up
    free(q);
    c_free_mechanism(linkage);
    for (i = 0; i < 4; ++i) c_free_mechanism_link(&links[i]);

    return 0;
}

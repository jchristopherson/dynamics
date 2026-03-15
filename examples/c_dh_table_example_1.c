#include <stdio.h>
#include <dynamics.h>

/*

This example builds the DH table for a 3R planar manipulator as shown in
Jazar's text as Example 126 on pg. 204.

*/
int main()
{
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
        R12[9], R123[9], T[16];

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
    ob3[0] = L1;    ob3[1] = 0.0;   ob3[2] = 0.0;

    c_matmul(3, 1, 3, 1.0, R1, 3, ob1, 3, 0.0, o1, 3);   // o1 = R1 * ob1
    c_matmul(3, 1, 3, 1.0, R12, 3, ob2, 3, 0.0, o2, 3);  // o2 = R1 * R2 * ob2
    c_matmul(3, 1, 3, 1.0, R123, 3, ob3, 3, 0.0, o3, 3); // o3 = R1 * R2 * R3 * ob3

    for (ii = 0; ii < 3; ++ii)
    {
        o2[ii] += o1[ii];
        o3[ii] += o2[ii];
    }

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

    // Display the contents of the table
    printf("Table Row Count: %i\n", tbl.count);
    printf("Frame\tLength\tTwist\tOffset\tAngle\n");
    printf("%i\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", 1,
        tbl.parameters[0].link_length,
        tbl.parameters[0].link_twist,
        tbl.parameters[0].link_offset,
        tbl.parameters[0].joint_angle
    );
    printf("%i\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", 2,
        tbl.parameters[1].link_length,
        tbl.parameters[1].link_twist,
        tbl.parameters[1].link_offset,
        tbl.parameters[1].joint_angle
    );
    printf("%i\t%0.3f\t%0.3f\t%0.3f\t%0.3f\n", 3,
        tbl.parameters[2].link_length,
        tbl.parameters[2].link_twist,
        tbl.parameters[2].link_offset,
        tbl.parameters[2].joint_angle
    );

    // Build the linkage transformation matrix
    c_dh_forward_kinematics_table(&tbl, T, 4);
    printf("\nLinkage Transformation Matrix:\n");
    printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\n", T[0], T[4], T[8], T[12]);
    printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\n", T[1], T[5], T[9], T[13]);
    printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\n", T[2], T[6], T[10], T[14]);
    printf("%0.3f\t%0.3f\t%0.3f\t%0.3f\n", T[3], T[7], T[11], T[15]);

    // Illustrate R123 for comparison
    printf("\nR123 Matrix:\n");
    printf("%0.3f\t%0.3f\t%0.3f\n", R123[0], R123[3], R123[6]);
    printf("%0.3f\t%0.3f\t%0.3f\n", R123[1], R123[4], R123[7]);
    printf("%0.3f\t%0.3f\t%0.3f\n", R123[2], R123[5], R123[8]);

    // Clean up
    c_free_dh_table(&tbl);

    // End
    return 0;
}
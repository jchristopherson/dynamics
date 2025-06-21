#include "dynamics_rotation_tests.h"
#include "dynamics_c_test_helper.h"
#include "dynamics.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

bool c_test_rotation_x()
{
    // Local Variables
    bool rst;
    double angle, r[9], ans[9];
    const double tol = 1.0e-8;
    
    // Initialization
    rst = true;
    angle = ((double)rand()) / RAND_MAX;
    ans[0] = 1.0;
    ans[1] = 0.0;
    ans[2] = 0.0;
    ans[3] = 0.0;
    ans[4] = cos(angle);
    ans[5] = sin(angle);
    ans[6] = 0.0;
    ans[7] = -sin(angle);
    ans[8] = cos(angle);

    // Test
    c_rotate_x(angle, r, 3);
    if (!compare_matrices(3, 3, r, 3, ans, 3, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_rotation_x -1");
    }

    // End
    return rst;
}


bool c_test_rotation_y()
{
    // Local Variables
    bool rst;
    double angle, r[9], ans[9];
    const double tol = 1.0e-8;
    
    // Initialization
    rst = true;
    angle = ((double)rand()) / RAND_MAX;
    ans[0] = cos(angle);
    ans[1] = 0.0;
    ans[2] = -sin(angle);
    ans[3] = 0.0;
    ans[4] = 1.0;
    ans[5] = 0.0;
    ans[6] = sin(angle);
    ans[7] = 0.0;
    ans[8] = cos(angle);

    // Test
    c_rotate_y(angle, r, 3);
    if (!compare_matrices(3, 3, r, 3, ans, 3, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_rotation_y -1");
    }

    // End
    return rst;
}


bool c_test_rotation_z()
{
    // Local Variables
    bool rst;
    double angle, r[9], ans[9];
    const double tol = 1.0e-8;
    
    // Initialization
    rst = true;
    angle = ((double)rand()) / RAND_MAX;
    ans[0] = cos(angle);
    ans[1] = sin(angle);
    ans[2] = 0.0;
    ans[3] = -sin(angle);
    ans[4] = cos(angle);
    ans[5] = 0.0;
    ans[6] = 0.0;
    ans[7] = 0.0;
    ans[8] = 1.0;

    // Test
    c_rotate_z(angle, r, 3);
    if (!compare_matrices(3, 3, r, 3, ans, 3, tol))
    {
        rst = false;
        printf("TEST FAILED: c_test_rotation_z -1");
    }

    // End
    return rst;
}

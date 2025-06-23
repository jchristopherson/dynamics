#include "dynamics.h"
#include "dynamics_c_test_helper.h"
#include "dynamics_stability_tests.h"
#include <stdio.h>
#include <math.h>

bool c_test_determine_local_stability()
{
    // Local Variables
    bool rst;
    const double d = 1.0e-1;
    const double a = 1.0;
    const double b = -5.0;
    double pt1, pt2, pt3, A1[4], A2[4], A3[4];
    int flag1, flag2, flag3;
    double complex ev1[2], ev2[2], ev3[2];

    // Initialization
    rst = true;
    pt1 = 0.0;
    pt2 = sqrt(-a / b);
    pt3 = -pt2;

    A1[0] = 0.0;
    A1[1] = -a - 3.0 * b * pt1 * pt1;
    A1[2] = 1.0;
    A1[3] = -d;

    A2[0] = 0.0;
    A2[1] = -a - 3.0 * b * pt2 * pt2;
    A2[2] = 1.0;
    A2[3] = -d;

    A3[0] = 0.0;
    A3[1] = -a - 3.0 * b * pt3 * pt3;
    A3[2] = 1.0;
    A3[3] = -d;

    // Test 1
    c_determine_local_stability(2, A1, 2, ev1, &flag1);
    if (flag1 != DYN_HYPERBOLIC_FIXED_POINT_SINK)
    {
        rst = false;
        printf("TEST FAILED: c_test_determine_local_stability -1\n");
    }

    // Test 2
    c_determine_local_stability(2, A2, 2, ev2, &flag2);
    if (flag2 != DYN_HYPERBOLIC_FIXED_POINT_SADDLE)
    {
        rst = false;
        printf("TEST FAILED: c_test_determine_local_stability -2\n");
    }

    // Test 3
    c_determine_local_stability(2, A3, 2, ev3, &flag3);
    if (flag3 != DYN_HYPERBOLIC_FIXED_POINT_SADDLE)
    {
        rst = false;
        printf("TEST FAILED: c_test_determine_local_stability -3\n");
    }

    // End
    return rst;
}

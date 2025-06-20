#include "dynamics_c_test_helper.h"
#include <math.h>
#include <stdbool.h>

bool compare_arrays(int n, const double *x, const double *y, double tol)
{
    // Local Variables
    bool rst;
    int i;

    // Process
    rst = true;
    for (i = 0; i < n; ++i)
    {
        if (fabs(x[i] - y[i]) > tol)
        {
            rst = false;
            break;
        }
    }

    // End
    return rst;
}
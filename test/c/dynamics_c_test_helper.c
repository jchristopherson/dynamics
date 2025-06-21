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

bool compare_matrices(int m, int n, const double *x, int ldx, const double *y,
    int ldy, double tol)
{
    // Local Variables
    int i, j;

    // Process
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            if (fabs(x[INDEX(i,j,ldx)] - y[INDEX(i,j,ldy)]) > tol)
            {
                return false;
            }
        }
    }
    return true;
}

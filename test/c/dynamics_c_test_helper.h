#ifndef DYNAMICS_C_TEST_HELPER_H_
#define DYNAMICS_C_TEST_HELPER_H_

#include <stdbool.h>

#define INDEX(i,j,m)((j) * (m) + (i))

bool compare_arrays(int n, const double *x, const double *y, double tol);
bool compare_matrices(int m, int n, const double *x, int ldx, const double *y,
    int ldy, double tol);
    
#endif

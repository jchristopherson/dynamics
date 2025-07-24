#ifndef DYNAMICS_C_TEST_HELPER_H_
#define DYNAMICS_C_TEST_HELPER_H_

#include <stdbool.h>

#define INDEX(i,j,m)((j) * (m) + (i))

bool compare_arrays(int n, const double *x, const double *y, double tol);
bool compare_matrices(int m, int n, const double *x, int ldx, const double *y,
    int ldy, double tol);
void mtx_mult(int m, int n, int p, const double *a, int lda, const double *b,
    int ldb, double *c, int ldc);
void transpose(int m, int n, const double *x, int ldx, double *xt, int ldxt);
void mtx_subtract(int m, int n, const double *x, int ldx, const double * y,
    int ldy, double *z, int ldz);
void create_random_vector(int n, double *x);
void to_skew_symmetric_matrix(const double *a, double *b, int ldb);
void zero_matrix(int m, int n, double *x, int ldx);
void vec_subtract(int n, const double *x, const double *y, double *z);
    
#endif

#include "dynamics_c_test_helper.h"
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>
#include <complex.h>

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

bool compare_complex_arrays(int n, const double complex *x, 
    const double complex *y, double tol)
{
    // Local Variables
    bool rst;
    int i;

    // Process
    rst = true;
    for (i = 0; i < n; ++i)
    {
        if (fabs(creal(x[i]) - creal(y[i])) > tol || 
            fabs(cimag(x[i]) - cimag(y[i])) > tol)
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

void mtx_mult(int m, int n, int p, const double *a, int lda, const double *b,
    int ldb, double *c, int ldc)
{
    // Local Variables
    int i, j, k;
    double val;

    // Process
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            val = 0.0;
            for (k = 0; k < p; ++k)
            {
                val += a[INDEX(i,k,lda)] * b[INDEX(k,j,ldb)];
            }
            c[INDEX(i,j,ldc)] = val;
        }
    }
}

void transpose(int m, int n, const double *x, int ldx, double *xt, int ldxt)
{
    // Local Variables
    int i, j;

    // Process
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            xt[INDEX(j,i,ldxt)] = x[INDEX(i,j,ldx)];
        }
    }
}

void mtx_subtract(int m, int n, const double *x, int ldx, const double * y,
    int ldy, double *z, int ldz)
{
    // Local Variables
    int i, j;

    // Process
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            z[INDEX(i,j,ldz)] = x[INDEX(i,j,ldx)] - y[INDEX(i,j,ldy)];
        }
    }
}


void create_random_vector(int n, double *x)
{
    // Local Variables
    int i;

    // Process
    for (i = 0; i < n; ++i)
    {
        x[i] = ((double)rand()) / RAND_MAX;
    }
}


void to_skew_symmetric_matrix(const double *a, double *b, int ldb)
{
    b[INDEX(0,0,ldb)] = 0.0;
    b[INDEX(1,0,ldb)] = a[2];
    b[INDEX(2,0,ldb)] = -a[1];
    b[INDEX(0,1,ldb)] = -a[2];
    b[INDEX(1,1,ldb)] = 0.0;
    b[INDEX(2,1,ldb)] = a[0];
    b[INDEX(0,2,ldb)] = a[1];
    b[INDEX(1,2,ldb)] = -a[0];
    b[INDEX(2,2,ldb)] = 0.0;
}


void zero_matrix(int m, int n, double *x, int ldx)
{
    // Local Variables
    int i, j;

    // Process
    for (j = 0; j < n; ++j)
    {
        for (i = 0; i < m; ++i)
        {
            x[INDEX(i,j,ldx)] = 0.0;
        }
    }
}

void vec_subtract(int n, const double *x, const double *y, double *z)
{
    // Local Variables
    int i;

    // Process
    for (i = 0; i < n; ++i)
    {
        z[i] = x[i] - y[i];
    }
}

bool compare_quaternion_to_array(const double x[4], const c_quaternion *q, double tol)
{
    if (fabs(x[0] - q->w) > tol) return false;
    if (fabs(x[1] - q->x) > tol) return false;
    if (fabs(x[2] - q->y) > tol) return false;
    if (fabs(x[3] - q->z) > tol) return false;
    return true;
}

bool compare_quaternions(const c_quaternion *x, const c_quaternion *y, double tol)
{
    if (fabs(x->w - y->w) > tol) return false;
    if (fabs(x->x - y->x) > tol) return false;
    if (fabs(x->y - y->y) > tol) return false;
    if (fabs(x->z - y->z) > tol) return false;
    return true;
}

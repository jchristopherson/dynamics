#include "dynamics.h"
#include <stdlib.h>


int alloc_dynamic_system_measurement(int n, c_dynamic_system_measurement *x);
void free_dynamic_system_measurement(c_dynamic_system_measurement *x);



int alloc_dynamic_system_measurement(int n, c_dynamic_system_measurement *x)
{
    size_t sz;
    sz = (size_t)(n * sizeof(double));
    x->npts = n;
    x->input = (double*)malloc(sz);
    if (!x->input) return -1;
    x->output = (double*)malloc(sz);
    if (!x->output) 
    {
        free_dynamic_system_measurement(x);
        return -1;
    }
    x->t = (double*)malloc(sz);
    if (!x->t)
    {
        free_dynamic_system_measurement(x);
        return -1;
    }
    return 0;
}

void free_dynamic_system_measurement(c_dynamic_system_measurement *x)
{
    x->npts = 0;
    if (x->input) free(x->input);
    if (x->output) free(x->output);
    if (x->t) free(x->t);
    x->input = NULL;
    x->output = NULL;
    x->t = NULL;
}

c_dynamic_system_measurement* alloc_dynamic_system_measurement_array(int n, 
    const int *ptsper)
{
    int i, flag;
    c_dynamic_system_measurement *obj;
    obj = (c_dynamic_system_measurement*)malloc(
        (size_t)(n * sizeof(c_dynamic_system_measurement))
    );
    for (i = 0; i < n; ++i)
    {
        flag = alloc_dynamic_system_measurement(ptsper[i], &obj[i]);
        if (flag != 0) 
        {
            free_dynamic_system_measurement_array(i, obj);
            return NULL;
        }
    }
    return obj;
}

void free_dynamic_system_measurement_array(int n, 
    c_dynamic_system_measurement *x)
{
    int i;
    if (!x) return;
    for (i = 0; i < n; ++i)
    {
        free_dynamic_system_measurement(&x[i]);
    }
    free(x);
    x = NULL;
}

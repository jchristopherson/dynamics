#include "dynamics.h"
#include <stdlib.h>


int c_alloc_dynamic_system_measurement(int n, c_dynamic_system_measurement *x)
{
    size_t sz;
    sz = (size_t)(n * sizeof(double));
    x->npts = n;
    x->input = (double*)malloc(sz);
    if (!x->input) return -1;
    x->output = (double*)malloc(sz);
    if (!x->output) 
    {
        c_free_dynamic_system_measurement(x);
        return -1;
    }
    x->t = (double*)malloc(sz);
    if (!x->t)
    {
        c_free_dynamic_system_measurement(x);
        return -1;
    }
    return 0;
}

void c_free_dynamic_system_measurement(c_dynamic_system_measurement *x)
{
    x->npts = 0;
    if (x->input) free(x->input);
    if (x->output) free(x->output);
    if (x->t) free(x->t);
    x->input = NULL;
    x->output = NULL;
    x->t = NULL;
}

c_dynamic_system_measurement* c_alloc_dynamic_system_measurement_array(int n, 
    const int *ptsper)
{
    int i, flag;
    c_dynamic_system_measurement *obj;
    obj = (c_dynamic_system_measurement*)malloc(
        (size_t)(n * sizeof(c_dynamic_system_measurement))
    );
    for (i = 0; i < n; ++i)
    {
        flag = c_alloc_dynamic_system_measurement(ptsper[i], &obj[i]);
        if (flag != 0) 
        {
            c_free_dynamic_system_measurement_array(n, obj);
            return NULL;
        }
    }
    return obj;
}

void c_free_dynamic_system_measurement_array(int n, 
    c_dynamic_system_measurement *x)
{
    int i;
    if (!x) return;
    for (i = 0; i < n; ++i)
    {
        c_free_dynamic_system_measurement(&x[i]);
    }
    free(x);
    x = NULL;
}

int c_alloc_dh_table(int n, c_dh_table *tbl)
{
    size_t sz;
    sz = (size_t)(n * sizeof(c_dh_parameter_set));
    tbl->count = n;
    tbl->parameters = (c_dh_parameter_set*)malloc(sz);
    return (!tbl->parameters) ? -1 : 0;
}

void c_free_dh_table(c_dh_table *tbl)
{
    tbl->count = 0;
    if (tbl->parameters) free(tbl->parameters);
    tbl->parameters = NULL;
}

int c_alloc_serial_linkage(int n, c_serial_linkage *lnk)
{
    size_t sz;
    sz = (size_t)(n * sizeof(c_binary_link));
    lnk->link_count = n;
    lnk->links = (c_binary_link*)malloc(sz);
    return (!lnk->links) ? -1 : 0;
}

void c_free_serial_linkage(c_serial_linkage *lnk)
{
    lnk->link_count = 0;
    if (lnk->links) free(lnk->links);
    lnk->links = NULL;
}

int c_alloc_polynomial(int order, c_polynomial *p)
{
    p->order = order;
    p->coefficients = (double*)malloc((size_t)(order + 1) * sizeof(double));
    return (!p->coefficients) ? -1 : 0;
}

void c_free_polynomial(c_polynomial *p)
{
    p->order = 0;
    if (p->coefficients) free(p->coefficients);
    p->coefficients = NULL;
}

int c_alloc_transfer_function(int numer_order, int denom_order, 
    c_transfer_function *tf)
{
    int flag;
    flag = c_alloc_polynomial(numer_order, &tf->numerator);
    if (flag < 0) return -1;
    flag = c_alloc_polynomial(denom_order, &tf->denominator);
    return flag;
}

void c_free_transfer_function(c_transfer_function *tf)
{
    c_free_polynomial(&tf->numerator);
    c_free_polynomial(&tf->denominator);
}

int c_alloc_state_space_model(int dimension, int n_inputs, int n_outputs, 
    c_state_space_model *mdl)
{
    mdl->dimension = dimension;
    mdl->n_inputs = n_inputs;
    mdl->n_outputs = n_outputs;
    mdl->A = (double*)malloc((size_t)(dimension * dimension * sizeof(double)));
    if (!mdl->A) return -1;
    mdl->B = (double*)malloc((size_t)(dimension * n_inputs * sizeof(double)));
    if (!mdl->B) 
    {
        free(mdl->A);
        mdl->A = NULL;
        return -1;
    }
    mdl->C = (double*)malloc((size_t)(n_outputs * dimension * sizeof(double)));
    if (!mdl->C)
    {
        free(mdl->A);
        mdl->A = NULL;
        free(mdl->B);
        mdl->B = NULL;
        return -1;
    }
    mdl->D = (double*)malloc((size_t)(n_outputs * n_inputs * sizeof(double)));
    if (!mdl->D)
    {
        free(mdl->A);
        mdl->A = NULL;
        free(mdl->B);
        mdl->B = NULL;
        free(mdl->C);
        mdl->C = NULL;
        return -1;
    }
    return 0;
}

void c_free_state_space_model(c_state_space_model *mdl)
{
    if (mdl->A) free(mdl->A);
    if (mdl->B) free(mdl->B);
    if (mdl->C) free(mdl->C);
    if (mdl->D) free(mdl->D);
    mdl->A = NULL;
    mdl->B = NULL;
    mdl->C = NULL;
    mdl->D = NULL;
    mdl->dimension = 0;
    mdl->n_inputs = 0;
    mdl->n_outputs = 0;
}

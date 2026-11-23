#include "dynamics.h"
#include <limits.h>
#include <stdint.h>
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

int c_alloc_mechanism_link(int nframes, c_mechanism_link *lnk)
{
    int i;
    size_t sz;
    sz = (size_t)(16 * nframes * sizeof(double));
    lnk->frame_count = nframes;
    lnk->frames = (double*)malloc(sz);
    if (!lnk->frames) return -1;
    for (i = 0; i < 16 * nframes; ++i) lnk->frames[i] = 0.0;
    for (i = 0; i < nframes; ++i)
    {
        lnk->frames[16*i] = 1.0;
        lnk->frames[16*i + 5] = 1.0;
        lnk->frames[16*i + 10] = 1.0;
        lnk->frames[16*i + 15] = 1.0;
    }
    lnk->mass = 1.0;
    for (i = 0; i < 3; ++i) lnk->cg[i] = 0.0;
    for (i = 0; i < 9; ++i) lnk->inertia[i] = 0.0;
    lnk->inertia[0] = 1.0;
    lnk->inertia[4] = 1.0;
    lnk->inertia[8] = 1.0;
    return 0;
}

void c_free_mechanism_link(c_mechanism_link *lnk)
{
    lnk->frame_count = 0;
    if (lnk->frames) free(lnk->frames);
    lnk->frames = NULL;
}

int c_alloc_polynomial(int order, c_polynomial *p)
{
    p->order = order;
    p->coefficients = NULL;
    if (order < 0 || order == INT_MAX ||
        (size_t)order + 1 > SIZE_MAX / sizeof(double)) return -1;
    p->coefficients = malloc(((size_t)order + 1) * sizeof(double));
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
    *tf = (c_transfer_function){0};
    if (c_alloc_polynomial(numer_order, &tf->numerator) < 0 ||
        c_alloc_polynomial(denom_order, &tf->denominator) < 0)
    {
        c_free_transfer_function(tf);
        return -1;
    }
    return 0;
}

void c_free_transfer_function(c_transfer_function *tf)
{
    c_free_polynomial(&tf->numerator);
    c_free_polynomial(&tf->denominator);
}

static int alloc_state_matrix(int rows, int cols, double **values)
{
    if (rows < 0 || cols < 0) return -1;
    if (rows == 0 || cols == 0)
    {
        /* Preserve a valid address for legacy Fortran bindings. */
        *values = malloc(sizeof(double));
        return *values ? 0 : -1;
    }
    if ((size_t)rows > SIZE_MAX / sizeof(double) / (size_t)cols) return -1;
    *values = malloc((size_t)rows * (size_t)cols * sizeof(double));
    return *values ? 0 : -1;
}

int c_alloc_state_space_model(int dimension, int n_inputs, int n_outputs,
    c_state_space_model *mdl)
{
    *mdl = (c_state_space_model){0};
    if (alloc_state_matrix(dimension, dimension, &mdl->A) < 0 ||
        alloc_state_matrix(dimension, n_inputs, &mdl->B) < 0 ||
        alloc_state_matrix(n_outputs, dimension, &mdl->C) < 0 ||
        alloc_state_matrix(n_outputs, n_inputs, &mdl->D) < 0)
    {
        c_free_state_space_model(mdl);
        return -1;
    }
    mdl->dimension = dimension;
    mdl->n_inputs = n_inputs;
    mdl->n_outputs = n_outputs;
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

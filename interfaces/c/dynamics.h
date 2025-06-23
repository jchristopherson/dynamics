#ifndef DYNAMICS_H_
#define DYNAMICS_H_

#include <complex.h>

#define DYN_HYPERBOLIC_FIXED_POINT_SINK 100
#define HYPERBOLIC_FIXED_POINT_SOURCE 101
#define DYN_HYPERBOLIC_FIXED_POINT_SADDLE 102
#define DYN_NONHYPERBOLIC_FIXED_POINT_UNSTABLE 103
#define DYN_NONHYPERBOLIC_FIXED_POINT_NEUTRALLY_STABLE 104
#define DYN_NONHYPERBOLIC_FIXED_POINT_CENTER 105

#ifdef __cplusplus
extern "C" {
#endif

double c_q_factor(double zeta);
double c_estimate_bandwidth(double fn, double zeta);
double c_logarithmic_decrement(double x1, double x2, int n);
double c_damping_from_log_decrement(double delta);
void c_find_free_response_properties(int n, const double *t, const double *x,
    double s, int np, double *delta, double *fn, double *x1, double *x2,
    double *t1, double *t2);
double c_rise_time(double wn, double zeta);
double c_find_settling_amplitude(int n, const double *x);
double c_damping_from_fractional_overshoot(int n, const double *x);
void c_evaluate_step_response(int n, double wn, double zeta, double xs,
    const double *t, double *x);

void c_rotate_x(double angle, double *r, int ldr);
void c_rotate_y(double angle, double *r, int ldr);
void c_rotate_z(double angle, double *r, int ldr);
void c_rotate(const double *i, const double *j, const double *k, double *r, 
    int ldr);
void c_acceleration_transform(const double *alpha, const double *omega,
    const double *a, const double *x, double *r, int ldr);
void c_velocity_transform(const double *omega, const double *v, 
    const double *x, double *r, int ldr);

void c_determine_local_stability(int n, const double *a, int lda,
    double complex *ev, int *flag);

#ifdef __cplusplus
}
#endif
#endif
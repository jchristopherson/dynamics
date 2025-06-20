#ifndef DYNAMICS_H_
#define DYNAMICS_H_

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

#ifdef __cplusplus
}
#endif
#endif
#ifndef DYNAMICS_H_
#define DYNAMICS_H_

#include <complex.h>
#include <stdbool.h>

#define DYN_HYPERBOLIC_FIXED_POINT_SINK 100
#define DYN_HYPERBOLIC_FIXED_POINT_SOURCE 101
#define DYN_HYPERBOLIC_FIXED_POINT_SADDLE 102
#define DYN_NONHYPERBOLIC_FIXED_POINT_UNSTABLE 103
#define DYN_NONHYPERBOLIC_FIXED_POINT_NEUTRALLY_STABLE 104
#define DYN_NONHYPERBOLIC_FIXED_POINT_CENTER 105

#define DYN_REVOLUTE_JOINT 0
#define DYN_PRISMATIC_JOINT 1

#define DYN_FRF_ACCELERANCE_MODEL 1
#define DYN_FRF_RECEPTANCE_MODEL 2

#define DYN_RUNGE_KUTTA_23 10
#define DYN_RUNGE_KUTTA_45 11
#define DYN_RUNGE_KUTTA_853 12
#define DYN_ROSENBROCK 13
#define DYN_BDF 14
#define DYN_ADAMS 15

#define DYN_ACCELERANCE_MODEL 1
#define DYN_RECEPTANCE_MODEL 2

#define DYN_H1 1
#define DYN_H2 2

typedef void (*c_vecfcn)(int nvar, int neqn, const double *x, double *f);
typedef void (*c_modal_excite)(int n, double freq, double complex *f);
typedef void (*c_harmonic_ode)(int n, double freq, double t, const double *x,
    double *dxdt);
typedef double (*c_window_function)(int n, int bin);
typedef void (*c_constraint_equations)(int n, int neqn, int nparam, 
    const double *xg, const double *fg, const double *xc, const double *p,
    double *fc);
typedef void (*c_ode)(int n, double t, const double *x, double *dxdt);

typedef struct {
    bool converge_on_chng;
    bool converge_on_fcn;
    bool converge_on_zero_diff;
    int fcn_count;
    int gradient_count;
    int iter_count;
    int jacobian_count;
} c_iteration_behavior;

typedef struct {
    int cycle_count;
    int transient_cycles;
    int points_per_cycle;
} c_frequency_sweep_controls;

typedef struct {
    double change_in_solution_tolerance;
    double gradient_tolerance;
    double iteration_improvement_tolerance;
    double residual_tolerance;
    int max_function_evaluations;
    int max_iteration_between_updates;
    int max_iteration_count;
} c_iteration_controls;

typedef struct {
    double confidence_interval;
    double probability;
    double standard_error;
    double t_statistic;
} c_regression_statistics;

typedef struct {
    int npts;
    double *input;
    double *output;
    double *t;
} c_dynamic_system_measurement;

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

void c_dh_forward_kinematics(int n, const double *alpha, const double *a,
    const double *theta, const double *d, double *T, int ldt);
void c_dh_forward_kinematics_2(const double *T1, int ldt1, const double *T2,
    int ldt2, double *T, int ldt);
void c_dh_forward_kinematics_3(const double *T1, int ldt1, const double *T2,
    int ldt2, const double *T3, int ldt3, double *T, int ldt);
void c_dh_forward_kinematics_4(const double *T1, int ldt1, const double *T2,
    int ldt2, const double *T3, int ldt3, const double *T4, double *T, int ldt);
void c_dh_forward_kinematics_5(const double *T1, int ldt1, const double *T2,
    int ldt2, const double *T3, int ldt3, const double *T4, int ldt4,
    const double *T5, int ldt5, double *T, int ldt);
void c_dh_forward_kinematics_6(const double *T1, int ldt1, const double *T2,
    int ldt2, const double *T3, int ldt3, const double *T4, int ldt4,
    const double *T5, int ldt5, const double *T6, int ldt6, double *T, int ldt);
void c_dh_forward_kinematics_7(const double *T1, int ldt1, const double *T2,
    int ldt2, const double *T3, int ldt3, const double *T4, int ldt4,
    const double *T5, int ldt5, const double *T6, int ldt6, const double *T7,
    int ldt7, double *T, int ldt);
void c_dh_forward_kinematics_8(const double *T1, int ldt1, const double *T2,
    int ldt2, const double *T3, int ldt3, const double *T4, int ldt4,
    const double *T5, int ldt5, const double *T6, int ldt6, const double *T7,
    int ldt7, const double *T8, int ldt8, double *T, int ldt);
void c_dh_jacobian(int n, const double *alpha, const double *a, 
    const double *theta, const double *d, const int *jtypes, double *jac,
    int ldjac);
void c_dh_matrix(double alpha, double a, double theta, double d, double *T,
    int ldt);
void c_dh_rotate_x(double alpha, double *T, int ldt);
void c_dh_rotate_z(double theta, double *T, int ldt);
void c_dh_translate_x(double a, double *T, int ldt);
void c_dh_translate_z(double d, double *T, int ldt);
void c_jacobian_generating_vector(const double *d, const double *k, 
    const double *R, int ldr, int jtype, double *jvec);
void c_solve_inverse_kinematics(int njoints, int neqn, const c_vecfcn mdl,
    const double *qo, const double *constraints, double *jvar, double *resid,
    c_iteration_behavior *ib);

void c_frequency_response(int n, int nfreq, const double *mass, int ldm,
    const double *stiff, int ldk, double alpha, double beta, const double *freq,
    const c_modal_excite frc, double *modes, double *modeshapes, int ldms,
    double complex *rsp, int ldr);
double c_compute_modal_damping(double lambda, double alpha, double beta);
double c_chirp(double t, double amp, double span, double f1Hz, double f2Hz);
void c_modal_response(int n, const double *mass, int ldm, const double *stiff,
    int ldk, double *freqs, double *modeshapes, int ldms);
void c_normalize_mode_shapes(int n, double *x, int ldx);
void c_frf_sweep(int n, int nfreq, c_harmonic_ode fcn, const double *freq,
    const double *iv, int solver, double complex *rsp, int ldr, 
    const c_frequency_sweep_controls *opts);
void c_set_frequency_sweep_defaults(c_frequency_sweep_controls *x);
void c_evaluate_accelerance_frf_model(int n, int norder, const double *mdl,
    const double *omega, double complex *h);
void c_evaluate_receptance_frf_model(int n, int norder, const double *mdl,
    const double *omega, double complex *h);
void c_set_iteration_controls_defaults(c_iteration_controls *x);
void c_fit_frf(int n, int norder, int method, const double *freq,
    const double complex *rsp, const double *maxp, const double *minp,
    const c_iteration_controls *controls, double *mdl, 
    c_regression_statistics *stats);
void c_siso_frequency_response(int n, int nf, const double *x, const double *y,
    double fs, int winsize, c_window_function winfun, int method, double *freq,
    double complex *rsp);

void c_cross_product(const double *x, const double *y, double *z);
void c_to_skew_symmetric(const double *x, double *y, int ldy);

void c_siso_model_fit_least_squares(int nsets, int nparams, int neqns, 
    const c_ode fcn, const c_dynamic_system_measurement *x, const double *ic,
    double *p, int integrator, int ind, const double *maxp, const double *minp,
    const c_iteration_controls *controls, int nconstraints, const double *xc,
    const double *yc, const c_constraint_equations constraints, int nweights,
    const double *weights, c_regression_statistics *stats, 
    c_iteration_behavior *info);

c_dynamic_system_measurement* alloc_dynamic_system_measurement_array(int n, 
    const int *ptsper);
void free_dynamic_system_measurement_array(int n, 
    c_dynamic_system_measurement *x);

#ifdef __cplusplus
}
#endif
#endif
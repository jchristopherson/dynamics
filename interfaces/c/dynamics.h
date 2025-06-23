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

typedef void (*c_vecfcn)(int nvar, int neqn, const double *x, double *f);
typedef void (*c_modal_excite)(int n, double freq, double *f);

typedef struct {
    bool converge_on_chng;
    bool converge_on_fcn;
    bool converge_on_zero_diff;
    int fcn_count;
    int gradient_count;
    int iter_count;
    int jacobian_count;
} c_iteration_behavior;

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

#ifdef __cplusplus
}
#endif
#endif
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

#define DYN_LEVENBERG_MARQUARDT_UPDATE 1
#define DYN_QUADRATIC_UPDATE 2
#define DYN_NIELSEN_UPDATE 3

typedef void (*c_vecfcn)(int nvar, int neqn, const double *x, double *f);
typedef void (*c_modal_excite)(int n, double freq, double complex *f);
typedef void (*c_harmonic_ode)(int n, double freq, double t, const double *x,
    double *dxdt);
typedef double (*c_window_function)(int n, int bin);
typedef void (*c_constraint_equations)(int n, int neqn, int nparam, 
    const double *xg, const double *fg, const double *xc, const double *p,
    double *fc);
typedef void (*c_ode)(int n, int nparam, const double *mdl, double t, 
    const double *x, double F, double *dxdt);

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

typedef struct 
{
    double damping_decrease_factor;
    double damping_increase_factor;
    double finite_difference_step_size;
    int method;
} c_lm_solver_options;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Estimates the Q-factor for a vibratory system.  
 * 
 * The Q-factor is computed as \f$ Q = \frac{1}{2 \zeta} \f$.
 * 
 * @param zeta The damping ratio.
 * @return The Q-factor.
 */
double c_q_factor(double zeta);
/**
 * @brief Estimates the bandwidth of the resonant mode of a vibratory system. 
 * The bandwidth is the width of the range of frequencies for which the energy 
 * is at least half its peak value and is computed as \f$ \Delta f = 
 * \frac{f_n}{Q} \f$.
 * 
 * @param fn The resonant frequency. The units are not important; however, the 
 * units of the output will be the same as the units of this parameter.
 * @param zeta The damping ratio.
 * @return The bandwidth.
 */
double c_estimate_bandwidth(double fn, double zeta);
/**
 * @brief Computes the logarithmic decrement given the value of two peaks in 
 * the time history of the free vibratory response of the system. The 
 * logarithmic decrement is calculated as follows.
 * 
 * \f$ \delta = \frac{1}{N} \ln \left( \frac{x(t)}{x(t + N T)} \right) =  
 * \frac{1}{N} \ln \left( \frac{x_1}{x_2} \right) \f$
 * 
 * @param x1 The amplitude of the first peak.
 * @param x2 The amplitude of the second peak that occurs N periods after the 
 * first.
 * @param n The number of periods of oscillation seperating the two peaks.
 * @return The logarithmic decrement.
 */
double c_logarithmic_decrement(double x1, double x2, int n);
/**
 * @brief Computes the damping ratio from the logarithmic decrement 
 * \f$.  The damping ratio is related to the logarithmic decrement by the
 * following relationship.
 * 
 * /f$ \zeta = \frac{\delta}{\sqrt{4 \pi^2 + \delta^2}} /f$
 * 
 * @param delta The logarithmic decrement.
 * @return The damping ratio.
 */
double c_damping_from_log_decrement(double delta);
/**
 * @brief Given a free-response time history, this routine attempts to find the 
 * logarithmic decrement and resonant frequency of a vibratory system. The 
 * logarithmic decrement is estimated by finding successive peaks by means of 
 * peak detection.
 * 
 * @param n The number of values in the response array.
 * @param t An N-element array containing the values in time.
 * @param x An N-element array containing the response sampled at the time 
 * points given in @p t.
 * @param s The sensitivity control of the peak detection algorithm.  A good
 * starting point is 0.1% of the peak-peak of the signal.
 * @param np The number of periods of oscillation seperating the two peaks
 * used to compute the logarithmic decrement. 
 * @param delta Returns the logarithmic decrement estimate. If sufficient peaks 
 * cannot be located, the routine returns NaN.
 * @param fn Returns the damped resonant frequency in units of Hz, assuming that
 * the time values are in seconds. If the time units are not in seconds, the 
 * units will be cycle/unit time with unit time being the units in which @p t 
 * is supplied. If sufficient peaks cannot be located, the routine returns NaN.
 * @param x1 Returns the amplitude of the first peak. If sufficient peaks cannot
 * be located, the routine returns NaN.
 * @param x2 Returns the amplitude of the second peak. If sufficient peaks 
 * cannot be located, the routine returns NaN.
 * @param t1 Returns the time at which the first peak was located. If sufficient
 * peaks cannot be located, the routine returns NaN.
 * @param t2 Returns the time at which the second peak was located. If 
 * sufficient peaks cannot be located, the routine returns NaN.
 */
void c_find_free_response_properties(int n, const double *t, const double *x,
    double s, int np, double *delta, double *fn, double *x1, double *x2,
    double *t1, double *t2);
/**
 * @brief Computes the rise time for an underdamped, second-order system. The 
 * rise time is the time it takes for the system response to go from 0% to 100% 
 * of its final value and is given by the following relationship.
 * 
 * /f$ t_r = \frac{1}{\omega_d} \left( \pi - 
 * \arctan \frac{\sqrt{1 - zeta^2}}{\zeta} \right) /f$
 * 
 * @param wn The resonant frequency of the system, in rad/s.
 * @param zeta The damping ratio of the system.  This value must be less than 1
 * as this relationship is only valid for an underdamped system.
 * @return The rise time, in units of seconds.
 */
double c_rise_time(double wn, double zeta);
/**
 * @brief Estimates the settling amplitude for a step response.
 * 
 * @param n The number of values in @p x.
 * @param x An N-element array containing the step response of the system.
 * @return The settling amplitude of the step response.
 */
double c_find_settling_amplitude(int n, const double *x);
/**
 * @brief Employs the method of fractional overshoot to estimate the damping 
 * ratio from the response of a system to a step input. This method is useful 
 * for cases where the damping ratio is between approximately 0.5 to 0.8. In 
 * such range, the logarithmic decrement approach becomes less precise.
 * 
 * The fractional overshoot method locates the amplitude of the first peak of 
 * oscillation (\f$x_p\f$) and the settling amplitude (\f$x_f\f$), and the 
 * estimates the damping ratio as follows.
 * 
 * \f$ s = \frac{x_p - x_f}{x_f} \f$
 * \f$ \zeta = \frac{1}{\sqrt{1 + \left( \frac{pi}{\ln{s}} \right)^2}} \f$
 * 
 * @param n The number of values in @p x.
 * @param x An N-element array containing the step response of the system.
 * @return The estimated damping ratio.
 */
double c_damping_from_fractional_overshoot(int n, const double *x);
/**
 * @brief Evaluates the response of an underdamped single-degree-of-freedom, 
 * linear system to a step function of amplitude \f$X_s\f$.  The step function 
 * response of an underdamped linear SDOF system is given as follows.
 * 
 * \f$ \ddot{x} + 2 \zeta \omega_n \dot{x} + \omega_n^2 x = \frac{F(t)}{m} \f$
 * \f$ \frac{x(t)}{X_s} = 1 - e^{-\zeta \omega_n t} \left( 
 * \frac{\zeta \omega_n}{\omega_d} \sin{\omega_d t} + \cos{\omega_d t}
 * \right) \f$
 * where,
 * \f$ \omega_d = \omega_n \sqrt{1 - \zeta^2} \f$
 * and
 * \f$ X_s = \frac{F}{k} \f$
 */
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
    const c_iteration_controls *controls, const c_lm_solver_options *opts,
    int nconstraints, const double *xc, const double *yc, 
    const c_constraint_equations constraints, int nweights,
    const double *weights, c_regression_statistics *stats, 
    c_iteration_behavior *info);
void c_set_lm_solver_options_defaults(c_lm_solver_options *x);
int alloc_dynamic_system_measurement(int n, c_dynamic_system_measurement *x);
void free_dynamic_system_measurement(c_dynamic_system_measurement *x);
c_dynamic_system_measurement* alloc_dynamic_system_measurement_array(int n, 
    const int *ptsper);
void free_dynamic_system_measurement_array(int n, 
    c_dynamic_system_measurement *x);

#ifdef __cplusplus
}
#endif
#endif
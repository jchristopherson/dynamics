#ifndef DYNAMICS_H_
#define DYNAMICS_H_

/**
 * @file dynamics.h
 * @brief Public C interface to the DYNAMICS library.
 *
 * Matrices use column-major storage to match the Fortran implementation.  A
 * matrix with `m` rows and `n` columns has leading dimension `ld` and stores
 * element `(i, j)` at `j * ld + i`, using zero-based C indices.  Link, joint,
 * and frame indices passed to the mechanism API are one-based.
 *
 * Objects returned by an allocation or creation routine remain owned by the
 * caller and must be released with the corresponding `c_free_*` routine.
 * Unless documented otherwise, pointer arguments must refer to storage large
 * enough for the dimensions supplied to the routine.
 */

#include <complex.h>
#include <stdbool.h>

/** @defgroup dynamics_constants Public constants */
/**@{*/
#define DYN_HYPERBOLIC_FIXED_POINT_SINK 100
#define DYN_HYPERBOLIC_FIXED_POINT_SOURCE 101
#define DYN_HYPERBOLIC_FIXED_POINT_SADDLE 102
#define DYN_NONHYPERBOLIC_FIXED_POINT_UNSTABLE 103
#define DYN_NONHYPERBOLIC_FIXED_POINT_NEUTRALLY_STABLE 104
#define DYN_NONHYPERBOLIC_FIXED_POINT_CENTER 105

#define DYN_REVOLUTE_JOINT 0
#define DYN_PRISMATIC_JOINT 1
#define DYN_FIXED_JOINT 2
#define DYN_CYLINDRICAL_JOINT 3
#define DYN_UNIVERSAL_JOINT 4
#define DYN_SPHERICAL_JOINT 5

#define DYN_FRF_ACCELERANCE_MODEL 1
#define DYN_FRF_RECEPTANCE_MODEL 2

#define DYN_RUNGE_KUTTA_23 10
#define DYN_RUNGE_KUTTA_45 11
#define DYN_RUNGE_KUTTA_853 12
#define DYN_ROSENBROCK 13
#define DYN_BDF 14
#define DYN_ADAMS 15
#define DYN_KENNEDY_CARPENTER_4 16
#define DYN_KENNEDY_CARPENTER_5 17
#define DYN_TSITOURAS_5 18

#define DYN_ACCELERANCE_MODEL 1
#define DYN_RECEPTANCE_MODEL 2

#define DYN_H1 1
#define DYN_H2 2

#define DYN_LEVENBERG_MARQUARDT_UPDATE 1
#define DYN_QUADRATIC_UPDATE 2
#define DYN_NIELSEN_UPDATE 3

#define DYN_POINCARE_TWO_SIDED 0
#define DYN_POINCARE_ONE_SIDED_FROM_FRONT 1
#define DYN_POINCARE_ONE_SIDED_FROM_BACK 2
/**@}*/

/** @brief Nonlinear vector function callback. */
typedef void (*c_vecfcn)(int nvar, int neqn, const double *x, double *f);
/** @brief Modal force callback used by frequency-response routines. */
typedef void (*c_modal_excite)(int n, double freq, double complex *f);
/** @brief Harmonic ordinary-differential-equation callback. */
typedef void (*c_harmonic_ode)(int n, double freq, double t, const double *x,
    double *dxdt);
/** @brief Window function callback used by SISO frequency analysis. */
typedef double (*c_window_function)(int n, int bin);
/** @brief Constraint callback used by least-squares system identification. */
typedef void (*c_constraint_equations)(int n, int neqn, int nparam, 
    const double *xg, const double *fg, const double *xc, const double *p,
    double *fc);
/** @brief ODE model callback used by system identification. */
typedef void (*c_ode_fit)(int n, int nparam, const double *mdl, double t, 
    const double *x, double F, double *dxdt);

/** @brief State-space input callback used by `c_lti_solve`. */
typedef void (*c_ss_excitation)(int n, double t, double *u);

/** @brief Iteration statistics returned by nonlinear solver routines. */
typedef struct {
    bool converge_on_chng;
    bool converge_on_fcn;
    bool converge_on_zero_diff;
    int fcn_count;
    int gradient_count;
    int iter_count;
    int jacobian_count;
} c_iteration_behavior;

/** @brief Controls for a frequency sweep. */
typedef struct {
    int cycle_count;
    int transient_cycles;
    int points_per_cycle;
    bool frequency_in_hz;
} c_frequency_sweep_controls;

/** @brief Controls for nonlinear iteration routines. */
typedef struct {
    double change_in_solution_tolerance;
    double gradient_tolerance;
    double iteration_improvement_tolerance;
    double residual_tolerance;
    int max_function_evaluations;
    int max_iteration_between_updates;
    int max_iteration_count;
} c_iteration_controls;

/** @brief Regression statistics returned by fitting routines. */
typedef struct {
    double confidence_interval;
    double probability;
    double standard_error;
    double t_statistic;
} c_regression_statistics;

/** @brief A single dynamic-system input/output measurement record. */
typedef struct {
    int npts;
    double *input;
    double *output;
    double *t;
} c_dynamic_system_measurement;

/** @brief Options for the Levenberg-Marquardt solver. */
typedef struct 
{
    double damping_decrease_factor;
    double damping_increase_factor;
    double finite_difference_step_size;
    int method;
} c_lm_solver_options;

/** @brief A quaternion stored as scalar component followed by vector terms. */
typedef struct
{
    double w;
    double x;
    double y;
    double z;
} c_quaternion;

/** @brief Plane coefficients satisfying `a*x + b*y + c*z + d = 0`. */
typedef struct
{
    double a;
    double b;
    double c;
    double d;
} c_plane;

/** @brief A line represented by a point and direction vector. */
typedef struct
{
    double r0[3];
    double v[3];
} c_line;

/** @brief A Pluecker line represented by direction and moment vectors. */
typedef struct
{
    double u[3];
    double m[3];
} c_plucker_line;

/** @brief An orthonormal coordinate system. */
typedef struct 
{
    double origin[3];
    double i[3];
    double j[3];
    double k[3];
} c_coordinate_system;

/** @brief One Denavit-Hartenberg parameter set. */
typedef struct
{
    double link_length;
    double link_twist;
    double link_offset;
    double joint_angle;
} c_dh_parameter_set;

/** @brief A dynamically allocated Denavit-Hartenberg table. */
typedef struct
{
    int count;
    c_dh_parameter_set *parameters;
} c_dh_table;

/** @brief One link in a serial Denavit-Hartenberg linkage. */
typedef struct
{
    double link_length;
    double link_twist;
    double link_offset;
    double joint_angle;
    int joint_type;     // DYN_REVOLUTE_JOINT or DYN_PRISMATIC_JOINT
    double mass;
    double cg[3];
    double inertia[9];  // 3-by-3 matrix in column-major format
} c_binary_link;

/** @brief A serial linkage and its link array. */
typedef struct
{
    int link_count;
    c_binary_link *links;
} c_serial_linkage;

/** @brief A multi-frame link used by a closed-loop mechanism. */
typedef struct
{
    int frame_count;
    double *frames;     // 4-by-4-by-frame_count array in column-major format
    double mass;
    double cg[3];
    double inertia[9];  // 3-by-3 matrix in column-major format
} c_mechanism_link;

/** @brief A joint connecting two mechanism link frames. */
typedef struct
{
    int joint_type;     // DYN_REVOLUTE_JOINT, DYN_PRISMATIC_JOINT, etc.
    int parent_link;    // one-based index of the parent link
    int parent_frame;   // one-based index of the frame on the parent link
    int child_link;     // one-based index of the child link
    int child_frame;    // one-based index of the frame on the child link
    bool actuated;
} c_joint;

/**
 * @brief Opaque handle to a closed-loop mechanism.
 *
 * Create handles with `c_create_parallel_linkage` or
 * `c_create_planar_linkage`, and release them with `c_free_mechanism`.
 */
typedef void* c_mechanism;

/** @brief A polynomial with dynamically allocated coefficients. */
typedef struct
{
    int order;
    double *coefficients;
} c_polynomial;

/** @brief A numerator/denominator transfer-function pair. */
typedef struct
{
    c_polynomial numerator;
    c_polynomial denominator;
} c_transfer_function;

/** @brief A continuous state-space model with column-major matrices. */
typedef struct
{
    int dimension;
    int n_inputs;
    int n_outputs;
    double *A;  // dimension -by- dimension
    double *B;  // dimension -by- n_inputs
    double *C;  // n_outputs -by- dimension
    double *D;  // n_outputs -by- n_inputs
} c_state_space_model;


#ifdef __cplusplus
extern "C" {
#endif

/** @defgroup dynamics_matrix Matrix and general kinematics */
/**@{*/
void c_matmul(int m, int n, int k, double alpha, const double *a, int lda,
    const double *b, int ldb, double beta, double *c, int ldc);

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
void c_rotate(const double i[3], const double j[3], const double k[3], double *r, 
    int ldr);
void c_acceleration_transform(const double alpha[3], const double omega[3],
    const double a[3], const double x[3], double *r, int ldr);
void c_velocity_transform(const double omega[3], const double v[3], 
    const double x[3], double *r, int ldr);

void c_determine_local_stability(int n, const double *a, int lda,
    double complex *ev, int *flag);

void c_dh_forward_kinematics_table(const c_dh_table *tbl, double *T, int ldt);
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
    const double *R, int ldr, int jtype, double jvec[6]);
void c_solve_inverse_kinematics(int njoints, int neqn, const c_vecfcn mdl,
    const double *qo, const double *constraints, const double *qmax,
    const double *qmin, double *jvar, double *resid,
    c_iteration_behavior *ib);
void c_to_angle_axis(const double *r, int ldr, double *angle, double axis[3]);
/**@}*/

/** @defgroup dynamics_frequency Frequency response and system identification */
/**@{*/
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

void c_cross_product(const double x[3], const double y[3], double z[3]);
void c_to_skew_symmetric(const double x[3], double *y, int ldy);
double c_vector_angle(const double x[3], const double y[3]);
double c_scalar_projection(const double x[3], const double y[3]);
void c_vector_projection(const double x[3], const double y[3], double z[3]);
double c_vector_magnitude(int n, const double *x);
void c_vector_normalize(int n, double *x);
double c_dot_product(int n, const double *x, const double *y);

void c_siso_model_fit_least_squares(int nsets, int nparams, int neqns, 
    const c_ode_fit fcn, const c_dynamic_system_measurement *x, 
    const double *ic, double *p, int integrator, int ind, const double *maxp, 
    const double *minp, const c_iteration_controls *controls, 
    const c_lm_solver_options *opts, int nconstraints, const double *xc, 
    const double *yc, const c_constraint_equations constraints, int nweights,
    const double *weights, c_regression_statistics *stats, 
    c_iteration_behavior *info);
void c_set_lm_solver_options_defaults(c_lm_solver_options *x);
int c_alloc_dynamic_system_measurement(int n, c_dynamic_system_measurement *x);
void c_free_dynamic_system_measurement(c_dynamic_system_measurement *x);
c_dynamic_system_measurement* c_alloc_dynamic_system_measurement_array(int n, 
    const int *ptsper);
void c_free_dynamic_system_measurement_array(int n, 
    c_dynamic_system_measurement *x);
/**@}*/

/** @defgroup dynamics_quaternion Quaternion operations */
/**@{*/
void c_quaternion_from_array(const double x[4], c_quaternion *q);
void c_quaternion_from_matrix(const double *x, int ldx, c_quaternion *q);
void c_quaternion_from_angle_axis(double angle, const double axis[3], 
    c_quaternion *q);
void c_quaternion_normalize(c_quaternion *q);
void c_quaternion_add(const c_quaternion *x, const c_quaternion *y,
    c_quaternion *q);
void c_quaternion_subtract(const c_quaternion *x, const c_quaternion *y,
    c_quaternion *q);
void c_quaternion_multiply(const c_quaternion *x, const c_quaternion *y,
    c_quaternion *q);   
void c_quaternion_divide(const c_quaternion *x, const c_quaternion *y,
    c_quaternion *q);   
void c_quaternion_scale(double x, const c_quaternion *y, c_quaternion *q);
void c_quaternion_conjugate(const c_quaternion *q, c_quaternion *qc);
void c_quaternion_rotate(const c_quaternion *q, const double r[3], double rp[3]);
double c_quaternion_abs(const c_quaternion *q);
void c_quaternion_inverse(const c_quaternion *q, c_quaternion *qinv);
void c_quaternion_to_matrix(const c_quaternion *q, double *r, int ldr);
void c_quaternion_to_angle_axis(const c_quaternion *q, double *angle, double axis[3]);
void c_quaternion_exp(const c_quaternion *q, c_quaternion *rst);
void c_quaternion_log(const c_quaternion *q, c_quaternion *rst);
void c_quaternion_pow(const c_quaternion *q, double exponent, c_quaternion *rst);
double c_quaternion_dot_product(const c_quaternion *x, const c_quaternion *y);
void c_quaternion_to_roll_pitch_yaw(const c_quaternion *q, double *roll, 
    double *pitch, double *yaw);
/**@}*/

/** @defgroup dynamics_geometry Geometry operations */
/**@{*/
void c_plane_normal(const c_plane* pln, double nrm[3]);
void c_plane_from_3_points(const double pt1[3], const double pt2[3], 
    const double pt3[3], c_plane *pln);
void c_plane_from_point_and_normal(const double pt[3], const double nrm[3],
    c_plane *pln);
void c_plane_from_points(int n, const double *pts, int ldp, c_plane *pln);
void c_flip_plane_normal(c_plane *pln);
void c_line_from_2_points(const double pt1[3], const double pt2[3], c_line *ln);
void c_line_from_2_planes(const c_plane *p1, const c_plane *p2, c_line *ln);
void c_line_from_points(int n, const double *pts, int ldp, c_line *ln);
void c_evaluate_line_position(const c_line *ln, double t, double x[3]);
bool c_is_parallel_vectors(int n, const double *x, const double *y, double tol);
bool c_is_parallel_lines(const c_line *x, const c_line *y, double tol);
bool c_is_parallel_planes(const c_plane *x, const c_plane *y, double tol);
bool c_is_point_on_plane(const double pt[3], const c_plane *pln, double tol);
bool c_is_point_on_line(const double pt[3], const c_line *ln, double tol);
double c_nearest_point_on_line(const double pt[3], const c_line *ln);
double c_point_to_line_distance(const double pt[3], const c_line *ln);
double c_point_to_plane_distance(const double pt[3], const c_plane *pln);
void c_vector_plane_projection(const double x[3], const c_plane *pln, double px[3]);
void c_point_plane_projection(const double pt[3], const c_plane *pln, double ppt[3]);
void c_plucker_line_from_2_points(const double pt1[3], const double pt2[3], 
    c_plucker_line *ln);
void c_plucker_line_from_line(const c_line *src, c_plucker_line *ln);
void c_plucker_line_from_2_planes(const c_plane *p1, const c_plane *p2, 
    c_plucker_line *ln);
void c_plucker_line_from_array(const double x[6], c_plucker_line *ln);
void c_plucker_line_mtx_mult(int n, const double *x, int ldx, 
    const c_plucker_line *ln, double *y);
void c_plucker_line_to_array(const c_plucker_line *ln, double x[6]);
void c_line_common_normal(const c_line *ln1, const c_line *ln2, c_line *ln);
void c_do_lines_intersect(const c_line *ln1, const c_line *ln2, bool *intersect,
    double *t1, double *t2, double tol);
void c_line_from_point_and_vector(const double pt[3], const double v[3],
    c_line *ln);

void c_poincare_map(int n, const double *x, const double *y, const double *z,
    const c_plane *pln, int side, int nbuffer, double *xbuffer, double *ybuffer,
    double *zbuffer, int *nactual);
/**@}*/

/** @defgroup dynamics_serial Serial linkage operations */
/**@{*/
int c_alloc_dh_table(int n, c_dh_table *tbl);
void c_free_dh_table(c_dh_table *tbl);
void c_define_link_csys(const double xim1[3], const double zim1[3], 
    const double zi[3], const double rim1[3], const double ri[3],
    c_coordinate_system *csys);
void c_define_csys(const double i[3], const double j[3], const double k[3], 
    const double o[3], c_coordinate_system *csys);
void c_build_dh_table(int n, const c_coordinate_system *csys, c_dh_table *tbl);

int c_alloc_serial_linkage(int n, c_serial_linkage *lnk);
void c_free_serial_linkage(c_serial_linkage *lnk);
void c_build_serial_linkage(int n, const c_binary_link *links, 
    c_serial_linkage *linkage);
void c_serial_linkage_forward_kinematics(int n, const c_serial_linkage *lnk,
    const double *q, double *T, int ldt);
void c_serial_linkage_jacobian(int n, const c_serial_linkage *lnk, 
    const double *q, double *jac, int ldj);
void c_serial_linkage_inverse_kinematics(int n, const c_serial_linkage *lnk,
    const double *qo, const double *trg, int ldt, double *q, 
    c_iteration_behavior *ib);
/**@}*/

/** @defgroup dynamics_parallel Parallel and planar linkage operations */
/**@{*/
int c_alloc_mechanism_link(int nframes, c_mechanism_link *lnk);
void c_free_mechanism_link(c_mechanism_link *lnk);
c_mechanism c_create_parallel_linkage(int nlinks, const c_mechanism_link *links,
    int njoints, const c_joint *joints, int base, int effector,
    const double *tool, int ldt);
c_mechanism c_create_planar_linkage(int nlinks, const c_mechanism_link *links,
    int njoints, const c_joint *joints, int base, int effector,
    const double *tool, int ldt);
void c_free_mechanism(c_mechanism obj);
int c_mechanism_link_count(c_mechanism obj);
int c_mechanism_joint_count(c_mechanism obj);
int c_mechanism_variable_count(c_mechanism obj);
int c_mechanism_loop_count(c_mechanism obj);
int c_mechanism_constraint_count(c_mechanism obj);
int c_mechanism_degrees_of_freedom(c_mechanism obj);
int c_mechanism_actuated_variable_count(c_mechanism obj);
int c_mechanism_space_dimension(c_mechanism obj);
int c_mechanism_link_frame_count(c_mechanism obj, int i);
void c_mechanism_link_frame(c_mechanism obj, int i, int k, double *T, int ldt);
void c_mechanism_get_configuration(c_mechanism obj, int n, double *q);
void c_mechanism_set_configuration(c_mechanism obj, int n, const double *q);
void c_mechanism_body_transform(c_mechanism obj, int i, int n, const double *q,
    double *T, int ldt);
void c_mechanism_end_effector_transform(c_mechanism obj, int n, 
    const double *q, double *T, int ldt);
void c_mechanism_constraints(c_mechanism obj, int n, const double *q, int nc,
    double *f);
void c_mechanism_constraint_jacobian(c_mechanism obj, int n, const double *q,
    double *jac, int ldj);
void c_mechanism_solve_configuration(c_mechanism obj, int na, const double *qa,
    int n, double *q, c_iteration_behavior *ib);
void c_mechanism_forward_kinematics(c_mechanism obj, int na, const double *qa,
    double *T, int ldt, c_iteration_behavior *ib);
void c_mechanism_jacobian(c_mechanism obj, int na, const double *qa, 
    double *jac, int ldj);
void c_mechanism_inverse_kinematics(c_mechanism obj, const double *trg, int ldt,
    int na, double *qa, c_iteration_behavior *ib);
/**@}*/

/** @defgroup dynamics_state Transfer functions and state-space models */
/**@{*/
int c_alloc_polynomial(int order, c_polynomial *p);
void c_free_polynomial(c_polynomial *p);
int c_alloc_transfer_function(int numer_order, int denom_order, 
    c_transfer_function *tf);
void c_free_transfer_function(c_transfer_function *tf);
int c_alloc_state_space_model(int dimension, int n_inputs, int n_outputs, 
    c_state_space_model *mdl);
void c_free_state_space_model(c_state_space_model *mdl);
void c_evaluate_transfer_function(const c_transfer_function *tf, int n,
    const double complex *s, double complex *z);
void c_transfer_function_poles(const c_transfer_function *tf, int n, 
    double complex *p);
void c_transfer_function_zeros(const c_transfer_function *tf, int n,
    double complex *z);
void c_to_ccf_state_space(const c_transfer_function *tf, c_state_space_model *ss);
void c_to_ocf_state_space(const c_transfer_function *tf, c_state_space_model *ss);
void c_create_state_space_model(int n, int n_out, const double *m, int ldm,
    const double *b, int ldb, const double *k, int ldk, 
    c_state_space_model *mdl);
void c_create_pid_state_space_model(double kp, double ki, double kd, double tau,
    const c_state_space_model *plant, c_state_space_model *mdl);
void c_transfer_function_multiply(const c_transfer_function *tf1,
    const c_transfer_function *tf2, c_transfer_function *tf);
void c_scale_transfer_function(double x, const c_transfer_function *tf1,
    c_transfer_function *tf);
void c_lti_solve(const c_state_space_model *mdl, const c_ss_excitation u,
    int n, const double *t, int ndof, const double *ic, int solver, 
    int nout, double *y, int ldy);
void c_state_space_poles(const c_state_space_model *mdl, int n,
    double complex *p);
void c_state_space_zeros(const c_state_space_model *mdl, int n,
    double complex *z, int *nz);
void c_state_space_transfer_function(const c_state_space_model *mdl, int nin,
    int nout, int n, const double complex *s, double complex *z, int ldz);
/**@}*/

#ifdef __cplusplus
}
#endif
#endif
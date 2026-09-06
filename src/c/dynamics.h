#ifndef DYNAMICS_H_
#define DYNAMICS_H_

/**
 * @file dynamics.h
 * @brief Public C interface to the DYNAMICS library. Matrices use column-major
 * storage to match the Fortran implementation. A matrix with `m` rows and `n`
 * columns has leading dimension `ld` and stores element `(i, j)` at `j * ld +
 * i`, using zero-based C indices. Link, joint, and frame indices passed to the
 * mechanism API are one-based. Objects returned by an allocation or creation
 * routine remain owned by the caller and must be released with the
 * corresponding `c_free_*` routine. Unless documented otherwise, pointer
 * arguments must refer to storage large enough for the dimensions supplied to
 * the routine.
 */

#include <complex.h>
#include <stdbool.h>

/**
 * @defgroup dynamics_constants Public constants
 */
/**
 * @{
 */
/**
 * Hyperbolic fixed point classified as a sink.
 */
#define DYN_HYPERBOLIC_FIXED_POINT_SINK 100
/**
 * Hyperbolic fixed point classified as a source.
 */
#define DYN_HYPERBOLIC_FIXED_POINT_SOURCE 101
/**
 * Hyperbolic fixed point classified as a saddle.
 */
#define DYN_HYPERBOLIC_FIXED_POINT_SADDLE 102
/**
 * Nonhyperbolic fixed point classified as unstable.
 */
#define DYN_NONHYPERBOLIC_FIXED_POINT_UNSTABLE 103
/**
 * Nonhyperbolic fixed point classified as neutrally stable.
 */
#define DYN_NONHYPERBOLIC_FIXED_POINT_NEUTRALLY_STABLE 104
/**
 * Nonhyperbolic fixed point classified as a center.
 */
#define DYN_NONHYPERBOLIC_FIXED_POINT_CENTER 105

/**
 * Revolute joint type.
 */
#define DYN_REVOLUTE_JOINT 0
/**
 * Prismatic joint type.
 */
#define DYN_PRISMATIC_JOINT 1
/**
 * Fixed joint type.
 */
#define DYN_FIXED_JOINT 2
/**
 * Cylindrical joint type.
 */
#define DYN_CYLINDRICAL_JOINT 3
/**
 * Universal joint type.
 */
#define DYN_UNIVERSAL_JOINT 4
/**
 * Spherical joint type.
 */
#define DYN_SPHERICAL_JOINT 5

/**
 * Frequency-response model expressed as accelerance.
 */
#define DYN_FRF_ACCELERANCE_MODEL 1
/**
 * Frequency-response model expressed as receptance.
 */
#define DYN_FRF_RECEPTANCE_MODEL 2

/**
 * Runge-Kutta 2/3 integration method.
 */
#define DYN_RUNGE_KUTTA_23 10
/**
 * Runge-Kutta 4/5 integration method.
 */
#define DYN_RUNGE_KUTTA_45 11
/**
 * Runge-Kutta 8/5/3 integration method.
 */
#define DYN_RUNGE_KUTTA_853 12
/**
 * Rosenbrock integration method.
 */
#define DYN_ROSENBROCK 13
/**
 * Backward-differentiation-formula integration method.
 */
#define DYN_BDF 14
/**
 * Adams integration method.
 */
#define DYN_ADAMS 15
/**
 * Kennedy-Carpenter fourth-order integration method.
 */
#define DYN_KENNEDY_CARPENTER_4 16
/**
 * Kennedy-Carpenter fifth-order integration method.
 */
#define DYN_KENNEDY_CARPENTER_5 17
/**
 * Tsitouras fifth-order integration method.
 */
#define DYN_TSITOURAS_5 18

/**
 * Accelerance model selector for system identification.
 */
#define DYN_ACCELERANCE_MODEL 1
/**
 * Receptance model selector for system identification.
 */
#define DYN_RECEPTANCE_MODEL 2

/**
 * H1 estimator selector for SISO frequency analysis.
 */
#define DYN_H1 1
/**
 * H2 estimator selector for SISO frequency analysis.
 */
#define DYN_H2 2

/**
 * Levenberg-Marquardt damping update method.
 */
#define DYN_LEVENBERG_MARQUARDT_UPDATE 1
/**
 * Quadratic damping update method.
 */
#define DYN_QUADRATIC_UPDATE 2
/**
 * Nielsen damping update method.
 */
#define DYN_NIELSEN_UPDATE 3

/**
 * Store intersections with both directions of the Poincare plane.
 */
#define DYN_POINCARE_TWO_SIDED 0
/**
 * Store intersections approaching the Poincare plane from the front.
 */
#define DYN_POINCARE_ONE_SIDED_FROM_FRONT 1
/**
 * Store intersections approaching the Poincare plane from the back.
 */
#define DYN_POINCARE_ONE_SIDED_FROM_BACK 2
/**
 * @}
 */

/**
 * Nonlinear vector function callback.
 * @param nvar Number of variables.
 * @param neqn Number of equations.
 * @param x Input variables.
 * @param f Output residual vector.
 */
typedef void (*c_vecfcn)(int nvar, int neqn, const double *x, double *f);
/**
 * Modal force callback used by frequency-response routines.
 * @param n Modal count.
 * @param freq Frequency.
 * @param f Output modal force.
 */
typedef void (*c_modal_excite)(int n, double freq, double complex *f);
/**
 * Harmonic ordinary-differential-equation callback.
 * @param n State dimension.
 * @param freq Excitation frequency.
 * @param t Time.
 * @param x State vector.
 * @param dxdt Output state derivative.
 */
typedef void (*c_harmonic_ode)(int n, double freq, double t, const double *x,
    double *dxdt);
/**
 * Window function callback used by SISO frequency analysis.
 * @param n Window length.
 * @param bin Zero-based sample index.
 * @return Window coefficient.
 */
typedef double (*c_window_function)(int n, int bin);
/**
 * Constraint callback used by least-squares system identification.
 * @param n Data-set index or count.
 * @param neqn Equation count.
 * @param nparam Parameter count.
 * @param xg General input values.
 * @param fg General constraint values.
 * @param xc Constraint inputs.
 * @param p Model parameters.
 * @param fc Output constraint residuals.
 */
typedef void (*c_constraint_equations)(int n, int neqn, int nparam, 
    const double *xg, const double *fg, const double *xc, const double *p,
    double *fc);
/**
 * ODE model callback used by system identification.
 * @param n State dimension.
 * @param nparam Parameter count.
 * @param mdl Model parameters.
 * @param t Time.
 * @param x State vector.
 * @param F Input or forcing value.
 * @param dxdt Output state derivative.
 */
typedef void (*c_ode_fit)(int n, int nparam, const double *mdl, double t, 
    const double *x, double F, double *dxdt);

/**
 * State-space input callback used by c_lti_solve.
 * @param n Input count.
 * @param t Time.
 * @param u Output input vector.
 */
typedef void (*c_ss_excitation)(int n, double t, double *u);

/**
 * @brief Iteration statistics returned by nonlinear solver routines.
 */
typedef struct {
    /**
     * True when convergence was detected from solution changes.
    */
    bool converge_on_chng;
    /**
     * True when convergence was detected from the function value.
    */
    bool converge_on_fcn;
    /**
     * True when convergence was detected from a zero difference.
    */
    bool converge_on_zero_diff;
    /**
     * Number of function evaluations.
    */
    int fcn_count;
    /**
     * Number of gradient evaluations.
    */
    int gradient_count;
    /**
     * Number of nonlinear iterations.
    */
    int iter_count;
    /**
     * Number of Jacobian evaluations.
    */
    int jacobian_count;
} c_iteration_behavior;

/**
 * @brief Controls for a frequency sweep.
 */
typedef struct {
    /**
     * Number of cycles in the frequency sweep.
    */
    int cycle_count;
    /**
     * Number of initial cycles discarded as transients.
    */
    int transient_cycles;
    /**
     * Number of frequency points evaluated per cycle.
    */
    int points_per_cycle;
    /**
     * True when sweep frequencies are specified in hertz.
    */
    bool frequency_in_hz;
} c_frequency_sweep_controls;

/**
 * @brief Controls for nonlinear iteration routines.
 */
typedef struct {
    /**
     * Convergence tolerance for changes in the solution.
    */
    double change_in_solution_tolerance;
    /**
     * Convergence tolerance for the gradient.
    */
    double gradient_tolerance;
    /**
     * Minimum improvement required between iterations.
    */
    double iteration_improvement_tolerance;
    /**
     * Convergence tolerance for the residual.
    */
    double residual_tolerance;
    /**
     * Maximum number of function evaluations.
    */
    int max_function_evaluations;
    /**
     * Maximum iterations between Jacobian or gradient updates.
    */
    int max_iteration_between_updates;
    /**
     * Maximum total number of iterations.
    */
    int max_iteration_count;
} c_iteration_controls;

/**
 * @brief Regression statistics returned by fitting routines.
 */
typedef struct {
    /**
     * Width of the confidence interval for a fitted parameter.
    */
    double confidence_interval;
    /**
     * Probability associated with the fitted parameter statistic.
    */
    double probability;
    /**
     * Standard error of the fitted parameter.
    */
    double standard_error;
    /**
     * Student t statistic for the fitted parameter.
    */
    double t_statistic;
} c_regression_statistics;

/**
 * @brief A single dynamic-system input/output measurement record.
 */
typedef struct {
    /**
     * Number of samples in the record.
    */
    int npts;
    /**
     * Input samples.
    */
    double *input;
    /**
     * Output samples.
    */
    double *output;
    /**
     * Sample times.
    */
    double *t;
} c_dynamic_system_measurement;

/**
 * @brief Options for the Levenberg-Marquardt solver.
 */
typedef struct 
{
    /**
     * Factor by which damping is decreased after a successful step.
    */
    double damping_decrease_factor;
    /**
     * Factor by which damping is increased after an unsuccessful step.
    */
    double damping_increase_factor;
    /**
     * Relative step used for finite-difference derivatives.
    */
    double finite_difference_step_size;
    /**
     * Damping update method, such as DYN_LEVENBERG_MARQUARDT_UPDATE.
    */
    int method;
} c_lm_solver_options;

/**
 * @brief A quaternion stored as scalar component followed by vector terms.
 */
typedef struct
{
    /**
     * Scalar quaternion component.
    */
    double w;
    /**
     * First vector quaternion component.
    */
    double x;
    /**
     * Second vector quaternion component.
    */
    double y;
    /**
     * Third vector quaternion component.
    */
    double z;
} c_quaternion;

/**
 * @brief Plane coefficients satisfying `a*x + b*y + c*z + d = 0`.
 */
typedef struct
{
    /**
     * Coefficient of x in the plane equation.
    */
    double a;
    /**
     * Coefficient of y in the plane equation.
    */
    double b;
    /**
     * Coefficient of z in the plane equation.
    */
    double c;
    /**
     * Constant plane coefficient.
    */
    double d;
} c_plane;

/**
 * @brief A line represented by a point and direction vector.
 */
typedef struct
{
    /**
     * A point on the line.
    */
    double r0[3];
    /**
     * Line direction vector.
    */
    double v[3];
} c_line;

/**
 * @brief A Pluecker line represented by direction and moment vectors.
 */
typedef struct
{
    /**
     * Pluecker line direction coordinates.
    */
    double u[3];
    /**
     * Pluecker line moment coordinates.
    */
    double m[3];
} c_plucker_line;

/**
 * @brief An orthonormal coordinate system.
 */
typedef struct 
{
    /**
     * Coordinate-system origin.
    */
    double origin[3];
    /**
     * First unit basis vector.
    */
    double i[3];
    /**
     * Second unit basis vector.
    */
    double j[3];
    /**
     * Third unit basis vector.
    */
    double k[3];
} c_coordinate_system;

/**
 * @brief One Denavit-Hartenberg parameter set.
 */
typedef struct
{
    /**
     * Link length parameter.
    */
    double link_length;
    /**
     * Link twist parameter.
    */
    double link_twist;
    /**
     * Link offset parameter.
    */
    double link_offset;
    /**
     * Joint angle parameter.
    */
    double joint_angle;
} c_dh_parameter_set;

/**
 * @brief A dynamically allocated Denavit-Hartenberg table.
 */
typedef struct
{
    /**
     * Number of parameter sets in the table.
    */
    int count;
    /**
     * Dynamically allocated parameter-set array.
    */
    c_dh_parameter_set *parameters;
} c_dh_table;

/**
 * @brief One link in a serial Denavit-Hartenberg linkage.
 */
typedef struct
{
    /**
     * Link length.
    */
    double link_length;
    /**
     * Link twist.
    */
    double link_twist;
    /**
     * Link offset.
    */
    double link_offset;
    /**
     * Joint angle or displacement.
    */
    double joint_angle;
    /**
     * Joint type, for example DYN_REVOLUTE_JOINT.
    */
    int joint_type;     // DYN_REVOLUTE_JOINT or DYN_PRISMATIC_JOINT
    /**
     * Link mass.
    */
    double mass;
    /**
     * Center-of-gravity position in the link frame.
    */
    double cg[3];
    /**
     * 3-by-3 inertia matrix in column-major storage.
    */
    double inertia[9];  // 3-by-3 matrix in column-major format
} c_binary_link;

/**
 * @brief A serial linkage and its link array.
 */
typedef struct
{
    /**
     * Number of links.
    */
    int link_count;
    /**
     * Link array.
    */
    c_binary_link *links;
} c_serial_linkage;

/**
 * @brief A multi-frame link used by a closed-loop mechanism.
 */
typedef struct
{
    /**
     * Number of frames attached to the link.
    */
    int frame_count;
    /**
     * 4-by-4-by-frame_count frame transforms in column-major storage.
    */
    double *frames;     // 4-by-4-by-frame_count array in column-major format
    /**
     * Link mass.
    */
    double mass;
    /**
     * Center-of-gravity position in the link frame.
    */
    double cg[3];
    /**
     * 3-by-3 inertia matrix in column-major storage.
    */
    double inertia[9];  // 3-by-3 matrix in column-major format
} c_mechanism_link;

/**
 * @brief A joint connecting two mechanism link frames.
 */
typedef struct
{
    /**
     * Joint type, such as DYN_REVOLUTE_JOINT or DYN_FIXED_JOINT.
    */
    int joint_type;     // DYN_REVOLUTE_JOINT, DYN_PRISMATIC_JOINT, etc.
    /**
     * One-based parent link index.
    */
    int parent_link;    // one-based index of the parent link
    /**
     * One-based parent-frame index.
    */
    int parent_frame;   // one-based index of the frame on the parent link
    /**
     * One-based child link index.
    */
    int child_link;     // one-based index of the child link
    /**
     * One-based child-frame index.
    */
    int child_frame;    // one-based index of the frame on the child link
    /**
     * Whether the joint variable is actuated.
    */
    bool actuated;
} c_joint;

/**
 * @brief Opaque handle to a closed-loop mechanism. Create handles with
 * `c_create_parallel_linkage` or `c_create_planar_linkage`, and release them
 * with `c_free_mechanism`.
 */
typedef void* c_mechanism;

/**
 * @brief A polynomial with dynamically allocated coefficients.
 */
typedef struct
{
    /**
     * Polynomial order.
    */
    int order;
    /**
     * Coefficients in ascending order of power.
    */
    double *coefficients;
} c_polynomial;

/**
 * @brief A numerator/denominator transfer-function pair.
 */
typedef struct
{
    /**
     * Numerator polynomial.
    */
    c_polynomial numerator;
    /**
     * Denominator polynomial.
    */
    c_polynomial denominator;
} c_transfer_function;

/**
 * @brief A continuous state-space model with column-major matrices.
 */
typedef struct
{
    /**
     * Number of state variables.
    */
    int dimension;
    /**
     * Number of inputs.
    */
    int n_inputs;
    /**
     * Number of outputs.
    */
    int n_outputs;
    /**
     * State matrix, dimension by dimension.
    */
    double *A;  // dimension -by- dimension
    /**
     * Input matrix, dimension by n_inputs.
    */
    double *B;  // dimension -by- n_inputs
    /**
     * Output matrix, n_outputs by dimension.
    */
    double *C;  // n_outputs -by- dimension
    /**
     * Feedthrough matrix, n_outputs by n_inputs.
    */
    double *D;  // n_outputs -by- n_inputs
} c_state_space_model;


#ifdef __cplusplus
extern "C" {
#endif

/**
 * @defgroup dynamics_matrix Matrix and general kinematics
 */
/**
 * @{
 */
/**
 * Multiply two column-major matrices: c = alpha*a*b + beta*c.
 */
/**
 * @param m Rows of a and c.
 * @param n Columns of b and c.
 * @param k Inner dimension.
 * @param alpha Scale factor for a*b.
 * @param a First matrix.
 * @param lda Leading dimension of a.
 * @param b Second matrix.
 * @param ldb Leading dimension of b.
 * @param beta Scale factor for c.
 * @param c Matrix updated in place.
 * @param ldc Leading dimension of c.
 */
void c_matmul(int m, int n, int k, double alpha, const double *a, int lda,
    const double *b, int ldb, double beta, double *c, int ldc);

/**
 * Compute the Q factor from a damping ratio.
 * @param zeta Damping ratio.
 * @return Q factor.
 */
double c_q_factor(double zeta);
/**
 * Estimate the half-power bandwidth.
 * @param fn Natural frequency.
 * @param zeta Damping ratio.
 * @return Bandwidth.
 */
double c_estimate_bandwidth(double fn, double zeta);
/**
 * Compute logarithmic decrement from two peaks.
 * @param x1 First peak.
 * @param x2 Second peak.
 * @param n Peak separation in cycles.
 * @return Logarithmic decrement.
 */
double c_logarithmic_decrement(double x1, double x2, int n);
/**
 * Convert logarithmic decrement to damping ratio.
 * @param delta Logarithmic decrement.
 * @return Damping ratio.
 */
double c_damping_from_log_decrement(double delta);
/**
 * Extract damping and frequency properties from a free response.
 * @param n Sample count.
 * @param t Sample times.
 * @param x Response samples.
 * @param s Settling threshold.
 * @param np Number of periods to use.
 * @param delta Output logarithmic decrement.
 * @param fn Output natural frequency.
 * @param x1 Output first peak.
 * @param x2 Output second peak.
 * @param t1 Output first peak time.
 * @param t2 Output second peak time.
 */
void c_find_free_response_properties(int n, const double *t, const double *x,
    double s, int np, double *delta, double *fn, double *x1, double *x2,
    double *t1, double *t2);
/**
 * Compute the 10-to-90 percent rise time.
 * @param wn Undamped natural frequency.
 * @param zeta Damping ratio.
 * @return Rise time.
 */
double c_rise_time(double wn, double zeta);
/**
 * Find the settled response amplitude.
 * @param n Sample count.
 * @param x Response samples.
 * @return Settling amplitude.
 */
double c_find_settling_amplitude(int n, const double *x);
/**
 * Estimate damping from fractional overshoot.
 * @param n Sample count.
 * @param x Response samples.
 * @return Damping ratio.
 */
double c_damping_from_fractional_overshoot(int n, const double *x);
/**
 * Evaluate a second-order step response.
 * @param n Sample count.
 * @param wn Undamped natural frequency.
 * @param zeta Damping ratio.
 * @param xs Step amplitude.
 * @param t Sample times.
 * @param x Output response samples.
 */
void c_evaluate_step_response(int n, double wn, double zeta, double xs,
    const double *t, double *x);

/**
 * Build a rotation matrix about the x axis.
 * @param angle Rotation angle in radians.
 * @param r Output 3-by-3 matrix.
 * @param ldr Leading dimension of r.
 */
void c_rotate_x(double angle, double *r, int ldr);
/**
 * Build a rotation matrix about the y axis.
 * @param angle Rotation angle in radians.
 * @param r Output 3-by-3 matrix.
 * @param ldr Leading dimension of r.
 */
void c_rotate_y(double angle, double *r, int ldr);
/**
 * Build a rotation matrix about the z axis.
 * @param angle Rotation angle in radians.
 * @param r Output 3-by-3 matrix.
 * @param ldr Leading dimension of r.
 */
void c_rotate_z(double angle, double *r, int ldr);
/**
 * Build a rotation matrix from three basis vectors.
 * @param i First basis vector.
 * @param j Second basis vector.
 * @param k Third basis vector.
 * @param r Output matrix.
 * @param ldr Leading dimension of r.
 */
void c_rotate(const double i[3], const double j[3], const double k[3], double *r, 
    int ldr);
/**
 * Build an acceleration transformation matrix.
 * @param alpha Angular acceleration.
 * @param omega Angular velocity.
 * @param a Linear acceleration.
 * @param x Position vector.
 * @param r Output matrix.
 * @param ldr Leading dimension of r.
 */
void c_acceleration_transform(const double alpha[3], const double omega[3],
    const double a[3], const double x[3], double *r, int ldr);
/**
 * Build a velocity transformation matrix.
 * @param omega Angular velocity.
 * @param v Linear velocity.
 * @param x Position vector.
 * @param r Output matrix.
 * @param ldr Leading dimension of r.
 */
void c_velocity_transform(const double omega[3], const double v[3], 
    const double x[3], double *r, int ldr);

/**
 * Classify the eigenvalues of a linearized system for local stability.
 * @param n Matrix order.
 * @param a System matrix.
 * @param lda Leading dimension of a.
 * @param ev Output eigenvalues.
 * @param flag Output stability classification.
 */
void c_determine_local_stability(int n, const double *a, int lda,
    double complex *ev, int *flag);

/**
 * Compute forward kinematics from a Denavit-Hartenberg table.
 * @param tbl Parameter table.
 * @param T Output transform.
 * @param ldt Leading dimension of T.
 */
void c_dh_forward_kinematics_table(const c_dh_table *tbl, double *T, int ldt);
/**
 * Compute serial Denavit-Hartenberg forward kinematics.
 * @param n Number of joints.
 * @param alpha Link twists.
 * @param a Link lengths.
 * @param theta Joint angles.
 * @param d Link offsets.
 * @param T Output transform.
 * @param ldt Leading dimension of T.
 */
void c_dh_forward_kinematics(int n, const double *alpha, const double *a,
    const double *theta, const double *d, double *T, int ldt);
/**
 * Multiply two homogeneous transforms.
 * @param T1 First transform.
 * @param ldt1 Leading dimension of T1.
 * @param T2 Second transform.
 * @param ldt2 Leading dimension of T2.
 * @param T Output transform.
 * @param ldt Leading dimension of T.
 */
void c_dh_forward_kinematics_2(const double *T1, int ldt1, const double *T2,
    int ldt2, double *T, int ldt);
/**
 * Multiply three homogeneous transforms.
 * @param T1 First transform.
 * @param ldt1 Leading dimension of T1.
 * @param T2 Second transform.
 * @param ldt2 Leading dimension of T2.
 * @param T3 Third transform.
 * @param ldt3 Leading dimension of T3.
 * @param T Output transform.
 * @param ldt Leading dimension of T.
 */
void c_dh_forward_kinematics_3(const double *T1, int ldt1, const double *T2,
    int ldt2, const double *T3, int ldt3, double *T, int ldt);
/**
 * Multiply four homogeneous transforms. Arguments follow
 * c_dh_forward_kinematics_3 and add T4/ldt4.
 */
void c_dh_forward_kinematics_4(const double *T1, int ldt1, const double *T2,
    int ldt2, const double *T3, int ldt3, const double *T4, double *T, int ldt);
/**
 * Multiply five homogeneous transforms. Arguments follow
 * c_dh_forward_kinematics_3 and add T4/ldt4 and T5/ldt5.
 */
void c_dh_forward_kinematics_5(const double *T1, int ldt1, const double *T2,
    int ldt2, const double *T3, int ldt3, const double *T4, int ldt4,
    const double *T5, int ldt5, double *T, int ldt);
/**
 * Multiply six homogeneous transforms. Arguments follow
 * c_dh_forward_kinematics_3 and add T4/ldt4, T5/ldt5, and T6/ldt6.
 */
void c_dh_forward_kinematics_6(const double *T1, int ldt1, const double *T2,
    int ldt2, const double *T3, int ldt3, const double *T4, int ldt4,
    const double *T5, int ldt5, const double *T6, int ldt6, double *T, int ldt);
/**
 * Multiply seven homogeneous transforms. Arguments follow
 * c_dh_forward_kinematics_3 and add T4/ldt4 through T7/ldt7.
 */
void c_dh_forward_kinematics_7(const double *T1, int ldt1, const double *T2,
    int ldt2, const double *T3, int ldt3, const double *T4, int ldt4,
    const double *T5, int ldt5, const double *T6, int ldt6, const double *T7,
    int ldt7, double *T, int ldt);
/**
 * Multiply eight homogeneous transforms. Arguments follow
 * c_dh_forward_kinematics_3 and add T4/ldt4 through T8/ldt8.
 */
void c_dh_forward_kinematics_8(const double *T1, int ldt1, const double *T2,
    int ldt2, const double *T3, int ldt3, const double *T4, int ldt4,
    const double *T5, int ldt5, const double *T6, int ldt6, const double *T7,
    int ldt7, const double *T8, int ldt8, double *T, int ldt);
/**
 * Compute a Denavit-Hartenberg linkage Jacobian.
 * @param n Joint count.
 * @param alpha Link twists.
 * @param a Link lengths.
 * @param theta Joint angles.
 * @param d Link offsets.
 * @param jtypes Joint types.
 * @param jac Output Jacobian.
 * @param ldjac Leading dimension of jac.
 */
void c_dh_jacobian(int n, const double *alpha, const double *a, 
    const double *theta, const double *d, const int *jtypes, double *jac,
    int ldjac);
/**
 * Build one Denavit-Hartenberg homogeneous transform.
 * @param alpha Link twist.
 * @param a Link length.
 * @param theta Joint angle.
 * @param d Link offset.
 * @param T Output transform.
 * @param ldt Leading dimension of T.
 */
void c_dh_matrix(double alpha, double a, double theta, double d, double *T,
    int ldt);
/**
 * Build the Denavit-Hartenberg x-axis rotation matrix.
 * @param alpha Rotation angle.
 * @param T Output matrix.
 * @param ldt Leading dimension of T.
 */
void c_dh_rotate_x(double alpha, double *T, int ldt);
/**
 * Build the Denavit-Hartenberg z-axis rotation matrix.
 * @param theta Rotation angle.
 * @param T Output matrix.
 * @param ldt Leading dimension of T.
 */
void c_dh_rotate_z(double theta, double *T, int ldt);
/**
 * Build the Denavit-Hartenberg x translation matrix.
 * @param a Translation distance.
 * @param T Output matrix.
 * @param ldt Leading dimension of T.
 */
void c_dh_translate_x(double a, double *T, int ldt);
/**
 * Build the Denavit-Hartenberg z translation matrix.
 * @param d Translation distance.
 * @param T Output matrix.
 * @param ldt Leading dimension of T.
 */
void c_dh_translate_z(double d, double *T, int ldt);
/**
 * Generate one joint Jacobian column.
 * @param d Joint position vector.
 * @param k Joint axis.
 * @param R Rotation matrix.
 * @param ldr Leading dimension of R.
 * @param jtype Joint type.
 * @param jvec Output six-vector.
 */
void c_jacobian_generating_vector(const double *d, const double *k, 
    const double *R, int ldr, int jtype, double jvec[6]);
/**
 * Solve a nonlinear inverse-kinematics problem.
 * @param njoints Number of joint variables.
 * @param neqn Number of equations.
 * @param mdl Residual callback.
 * @param qo Initial joint variables.
 * @param constraints Constraint values.
 * @param qmax Upper bounds.
 * @param qmin Lower bounds.
 * @param jvar Output joint variables.
 * @param resid Output residual.
 * @param ib Output iteration statistics.
 */
void c_solve_inverse_kinematics(int njoints, int neqn, const c_vecfcn mdl,
    const double *qo, const double *constraints, const double *qmax,
    const double *qmin, double *jvar, double *resid,
    c_iteration_behavior *ib);
/**
 * Convert a rotation matrix to angle-axis form.
 * @param r Rotation matrix.
 * @param ldr Leading dimension of r.
 * @param angle Output angle in radians.
 * @param axis Output unit axis.
 */
void c_to_angle_axis(const double *r, int ldr, double *angle, double axis[3]);
/**
 * @}
 */

/**
 * @defgroup dynamics_frequency Frequency response and system identification
 */
/**
 * @{
 */
/**
 * Compute a modal frequency response for a second-order system.
 * @param n System order.
 * @param nfreq Frequency count.
 * @param mass Mass matrix.
 * @param ldm Leading dimension of mass.
 * @param stiff Stiffness matrix.
 * @param ldk Leading dimension of stiff.
 * @param alpha Mass-proportional damping coefficient.
 * @param beta Stiffness-proportional damping coefficient.
 * @param freq Frequencies.
 * @param frc Modal force callback.
 * @param modes Output modal frequencies.
 * @param modeshapes Output mode shapes.
 * @param ldms Leading dimension of modeshapes.
 * @param rsp Output complex response.
 * @param ldr Leading dimension of rsp.
 */
void c_frequency_response(int n, int nfreq, const double *mass, int ldm,
    const double *stiff, int ldk, double alpha, double beta, const double *freq,
    const c_modal_excite frc, double *modes, double *modeshapes, int ldms,
    double complex *rsp, int ldr);
/**
 * Compute modal damping from Rayleigh coefficients.
 * @param lambda Modal eigenvalue.
 * @param alpha Mass-proportional coefficient.
 * @param beta Stiffness-proportional coefficient.
 * @return Modal damping ratio.
 */
double c_compute_modal_damping(double lambda, double alpha, double beta);
/**
 * Evaluate a swept-frequency chirp signal.
 * @param t Time.
 * @param amp Amplitude.
 * @param span Sweep span.
 * @param f1Hz Initial frequency in hertz.
 * @param f2Hz Final frequency in hertz.
 * @return Signal value.
 */
double c_chirp(double t, double amp, double span, double f1Hz, double f2Hz);
/**
 * Compute modal frequencies and mass-normalized mode shapes.
 * @param n System order.
 * @param mass Mass matrix.
 * @param ldm Leading dimension of mass.
 * @param stiff Stiffness matrix.
 * @param ldk Leading dimension of stiff.
 * @param freqs Output modal frequencies.
 * @param modeshapes Output mode shapes.
 * @param ldms Leading dimension of modeshapes.
 */
void c_modal_response(int n, const double *mass, int ldm, const double *stiff,
    int ldk, double *freqs, double *modeshapes, int ldms);
/**
 * Normalize mode-shape columns.
 * @param n Number of rows or modes.
 * @param x Mode-shape matrix, updated in place.
 * @param ldx Leading dimension of x.
 */
void c_normalize_mode_shapes(int n, double *x, int ldx);
/**
 * Perform a nonlinear harmonic frequency sweep.
 * @param n State dimension.
 * @param nfreq Frequency count.
 * @param fcn Harmonic ODE callback.
 * @param freq Frequencies.
 * @param iv Initial values.
 * @param solver Integration method.
 * @param rsp Output complex response.
 * @param ldr Leading dimension of rsp.
 * @param opts Sweep controls.
 */
void c_frf_sweep(int n, int nfreq, c_harmonic_ode fcn, const double *freq,
    const double *iv, int solver, double complex *rsp, int ldr, 
    const c_frequency_sweep_controls *opts);
/**
 * Fill frequency-sweep controls with defaults.
 * @param x Controls updated in place.
 */
void c_set_frequency_sweep_defaults(c_frequency_sweep_controls *x);
/**
 * Evaluate an accelerance rational FRF model.
 * @param n Frequency count.
 * @param norder Model order.
 * @param mdl Model coefficients.
 * @param omega Angular frequencies.
 * @param h Output complex FRF values.
 */
void c_evaluate_accelerance_frf_model(int n, int norder, const double *mdl,
    const double *omega, double complex *h);
/**
 * Evaluate a receptance rational FRF model.
 * @param n Frequency count.
 * @param norder Model order.
 * @param mdl Model coefficients.
 * @param omega Angular frequencies.
 * @param h Output complex FRF values.
 */
void c_evaluate_receptance_frf_model(int n, int norder, const double *mdl,
    const double *omega, double complex *h);
/**
 * Fill nonlinear iteration controls with defaults.
 * @param x Controls updated in place.
 */
void c_set_iteration_controls_defaults(c_iteration_controls *x);
/**
 * Fit an FRF model to measured data.
 * @param n Data count.
 * @param norder Model order.
 * @param method Model or update method.
 * @param freq Frequencies.
 * @param rsp Measured complex response.
 * @param maxp Parameter upper bounds.
 * @param minp Parameter lower bounds.
 * @param controls Iteration controls.
 * @param mdl Output fitted parameters.
 * @param stats Output regression statistics.
 */
void c_fit_frf(int n, int norder, int method, const double *freq,
    const double complex *rsp, const double *maxp, const double *minp,
    const c_iteration_controls *controls, double *mdl, 
    c_regression_statistics *stats);
/**
 * Estimate a SISO frequency response from input and output records.
 * @param n Sample count.
 * @param nf Frequency-bin count.
 * @param x Input samples.
 * @param y Output samples.
 * @param fs Sampling frequency.
 * @param winsize Window size.
 * @param winfun Window callback.
 * @param method H1 or H2 estimator.
 * @param freq Output frequencies.
 * @param rsp Output complex response.
 */
void c_siso_frequency_response(int n, int nf, const double *x, const double *y,
    double fs, int winsize, c_window_function winfun, int method, double *freq,
    double complex *rsp);

/**
 * Compute the cross product of two three-vectors.
 * @param x First vector.
 * @param y Second vector.
 * @param z Output vector.
 */
void c_cross_product(const double x[3], const double y[3], double z[3]);
/**
 * Form a skew-symmetric matrix from a three-vector.
 * @param x Input vector.
 * @param y Output matrix.
 * @param ldy Leading dimension of y.
 */
void c_to_skew_symmetric(const double x[3], double *y, int ldy);
/**
 * Compute the angle between two vectors.
 * @param x First vector.
 * @param y Second vector.
 * @return Angle in radians.
 */
double c_vector_angle(const double x[3], const double y[3]);
/**
 * Compute the scalar projection of one vector onto another.
 * @param x Vector being projected.
 * @param y Reference vector.
 * @return Scalar projection.
 */
double c_scalar_projection(const double x[3], const double y[3]);
/**
 * Project one vector onto another.
 * @param x Vector being projected.
 * @param y Reference vector.
 * @param z Output projection.
 */
void c_vector_projection(const double x[3], const double y[3], double z[3]);
/**
 * Compute a vector magnitude.
 * @param n Vector length.
 * @param x Vector.
 * @return Euclidean magnitude.
 */
double c_vector_magnitude(int n, const double *x);
/**
 * Normalize a vector in place.
 * @param n Vector length.
 * @param x Vector updated in place.
 */
void c_vector_normalize(int n, double *x);
/**
 * Compute a dot product.
 * @param n Vector length.
 * @param x First vector.
 * @param y Second vector.
 * @return Dot product.
 */
double c_dot_product(int n, const double *x, const double *y);

/**
 * Fit an ODE model to dynamic-system measurements by constrained least squares.
 * @param nsets Number of measurement sets.
 * @param nparams Number of parameters.
 * @param neqns Number of model equations.
 * @param fcn ODE callback.
 * @param x Measurement records.
 * @param ic Initial conditions.
 * @param p Parameters, updated in place.
 * @param integrator ODE integration method.
 * @param ind Integration direction or index.
 * @param maxp Upper parameter bounds.
 * @param minp Lower parameter bounds.
 * @param controls Iteration controls.
 * @param opts Levenberg-Marquardt options.
 * @param nconstraints Constraint count.
 * @param xc Constraint inputs.
 * @param yc Constraint outputs.
 * @param constraints Constraint callback.
 * @param nweights Weight count.
 * @param weights Residual weights.
 * @param stats Output regression statistics.
 * @param info Output iteration statistics.
 */
void c_siso_model_fit_least_squares(int nsets, int nparams, int neqns, 
    const c_ode_fit fcn, const c_dynamic_system_measurement *x, 
    const double *ic, double *p, int integrator, int ind, const double *maxp, 
    const double *minp, const c_iteration_controls *controls, 
    const c_lm_solver_options *opts, int nconstraints, const double *xc, 
    const double *yc, const c_constraint_equations constraints, int nweights,
    const double *weights, c_regression_statistics *stats, 
    c_iteration_behavior *info);
/**
 * Fill Levenberg-Marquardt options with defaults.
 * @param x Options updated in place.
 */
void c_set_lm_solver_options_defaults(c_lm_solver_options *x);
/**
 * Allocate one dynamic-system measurement record.
 * @param n Sample count.
 * @param x Record initialized by the routine.
 * @return Zero on success; nonzero on allocation failure.
 */
int c_alloc_dynamic_system_measurement(int n, c_dynamic_system_measurement *x);
/**
 * Release one dynamic-system measurement record.
 * @param x Record to release.
 */
void c_free_dynamic_system_measurement(c_dynamic_system_measurement *x);
/**
 * Allocate an array of measurement records.
 * @param n Number of records.
 * @param ptsper Samples per record.
 * @return Allocated array, or NULL on failure.
 */
c_dynamic_system_measurement* c_alloc_dynamic_system_measurement_array(int n, 
    const int *ptsper);
/**
 * Release an array of measurement records.
 * @param n Number of records.
 * @param x Array to release.
 */
void c_free_dynamic_system_measurement_array(int n, 
    c_dynamic_system_measurement *x);
/**
 * @}
 */

/**
 * @defgroup dynamics_quaternion Quaternion operations
 */
/**
 * @{
 */
/**
 * Construct a quaternion from four scalar components.
 * @param x Four components in w,x,y,z order.
 * @param q Output quaternion.
 */
void c_quaternion_from_array(const double x[4], c_quaternion *q);
/**
 * Construct a quaternion from a rotation matrix.
 * @param x Rotation matrix.
 * @param ldx Leading dimension of x.
 * @param q Output quaternion.
 */
void c_quaternion_from_matrix(const double *x, int ldx, c_quaternion *q);
/**
 * Construct a quaternion from angle-axis data.
 * @param angle Rotation angle in radians.
 * @param axis Unit rotation axis.
 * @param q Output quaternion.
 */
void c_quaternion_from_angle_axis(double angle, const double axis[3], 
    c_quaternion *q);
/**
 * Normalize a quaternion in place.
 * @param q Quaternion updated in place.
 */
void c_quaternion_normalize(c_quaternion *q);
/**
 * Add two quaternions.
 * @param x First quaternion.
 * @param y Second quaternion.
 * @param q Output quaternion.
 */
void c_quaternion_add(const c_quaternion *x, const c_quaternion *y,
    c_quaternion *q);
/**
 * Subtract two quaternions.
 * @param x First quaternion.
 * @param y Second quaternion.
 * @param q Output quaternion.
 */
void c_quaternion_subtract(const c_quaternion *x, const c_quaternion *y,
    c_quaternion *q);
/**
 * Multiply two quaternions.
 * @param x First quaternion.
 * @param y Second quaternion.
 * @param q Output quaternion.
 */
void c_quaternion_multiply(const c_quaternion *x, const c_quaternion *y,
    c_quaternion *q);   
/**
 * Divide two quaternions.
 * @param x Numerator quaternion.
 * @param y Denominator quaternion.
 * @param q Output quaternion.
 */
void c_quaternion_divide(const c_quaternion *x, const c_quaternion *y,
    c_quaternion *q);   
/**
 * Scale a quaternion.
 * @param x Scale factor.
 * @param y Input quaternion.
 * @param q Output quaternion.
 */
void c_quaternion_scale(double x, const c_quaternion *y, c_quaternion *q);
/**
 * Compute a quaternion conjugate.
 * @param q Input quaternion.
 * @param qc Output conjugate.
 */
void c_quaternion_conjugate(const c_quaternion *q, c_quaternion *qc);
/**
 * Rotate a vector using a quaternion.
 * @param q Rotation quaternion.
 * @param r Input vector.
 * @param rp Output rotated vector.
 */
void c_quaternion_rotate(const c_quaternion *q, const double r[3], double rp[3]);
/**
 * Compute the quaternion norm.
 * @param q Quaternion.
 * @return Quaternion magnitude.
 */
double c_quaternion_abs(const c_quaternion *q);
/**
 * Compute a quaternion inverse.
 * @param q Input quaternion.
 * @param qinv Output inverse.
 */
void c_quaternion_inverse(const c_quaternion *q, c_quaternion *qinv);
/**
 * Convert a quaternion to a rotation matrix.
 * @param q Quaternion.
 * @param r Output matrix.
 * @param ldr Leading dimension of r.
 */
void c_quaternion_to_matrix(const c_quaternion *q, double *r, int ldr);
/**
 * Convert a quaternion to angle-axis data.
 * @param q Quaternion.
 * @param angle Output angle.
 * @param axis Output axis.
 */
void c_quaternion_to_angle_axis(const c_quaternion *q, double *angle, double axis[3]);
/**
 * Compute the quaternion exponential.
 * @param q Input quaternion.
 * @param rst Output quaternion.
 */
void c_quaternion_exp(const c_quaternion *q, c_quaternion *rst);
/**
 * Compute the quaternion logarithm.
 * @param q Input quaternion.
 * @param rst Output quaternion.
 */
void c_quaternion_log(const c_quaternion *q, c_quaternion *rst);
/**
 * Raise a quaternion to a real power.
 * @param q Input quaternion.
 * @param exponent Power.
 * @param rst Output quaternion.
 */
void c_quaternion_pow(const c_quaternion *q, double exponent, c_quaternion *rst);
/**
 * Compute the quaternion dot product.
 * @param x First quaternion.
 * @param y Second quaternion.
 * @return Dot product.
 */
double c_quaternion_dot_product(const c_quaternion *x, const c_quaternion *y);
/**
 * Convert a quaternion to roll, pitch, and yaw angles.
 * @param q Quaternion.
 * @param roll Output roll.
 * @param pitch Output pitch.
 * @param yaw Output yaw.
 */
void c_quaternion_to_roll_pitch_yaw(const c_quaternion *q, double *roll, 
    double *pitch, double *yaw);
/**
 * @}
 */

/**
 * @defgroup dynamics_geometry Geometry operations
 */
/**
 * @{
 */
/**
 * Extract a plane normal.
 * @param pln Plane.
 * @param nrm Output normal vector.
 */
void c_plane_normal(const c_plane* pln, double nrm[3]);
/**
 * Construct a plane through three points.
 * @param pt1 First point.
 * @param pt2 Second point.
 * @param pt3 Third point.
 * @param pln Output plane.
 */
void c_plane_from_3_points(const double pt1[3], const double pt2[3], 
    const double pt3[3], c_plane *pln);
/**
 * Construct a plane from a point and normal.
 * @param pt Point on plane.
 * @param nrm Plane normal.
 * @param pln Output plane.
 */
void c_plane_from_point_and_normal(const double pt[3], const double nrm[3],
    c_plane *pln);
/**
 * Fit a plane to point data.
 * @param n Point count.
 * @param pts Point matrix.
 * @param ldp Leading dimension of pts.
 * @param pln Output plane.
 */
void c_plane_from_points(int n, const double *pts, int ldp, c_plane *pln);
/**
 * Reverse a plane normal and its equation.
 * @param pln Plane updated in place.
 */
void c_flip_plane_normal(c_plane *pln);
/**
 * Construct a line through two points.
 * @param pt1 First point.
 * @param pt2 Second point.
 * @param ln Output line.
 */
void c_line_from_2_points(const double pt1[3], const double pt2[3], c_line *ln);
/**
 * Construct the intersection line of two planes.
 * @param p1 First plane.
 * @param p2 Second plane.
 * @param ln Output line.
 */
void c_line_from_2_planes(const c_plane *p1, const c_plane *p2, c_line *ln);
/**
 * Fit a line to point data.
 * @param n Point count.
 * @param pts Point matrix.
 * @param ldp Leading dimension of pts.
 * @param ln Output line.
 */
void c_line_from_points(int n, const double *pts, int ldp, c_line *ln);
/**
 * Evaluate a point on a line.
 * @param ln Line.
 * @param t Line parameter.
 * @param x Output position.
 */
void c_evaluate_line_position(const c_line *ln, double t, double x[3]);
/**
 * Test whether two vectors are parallel within a tolerance.
 * @param n Vector length.
 * @param x First vector.
 * @param y Second vector.
 * @param tol Angular or residual tolerance.
 * @return true when parallel.
 */
bool c_is_parallel_vectors(int n, const double *x, const double *y, double tol);
/**
 * Test whether two lines are parallel.
 * @param x First line.
 * @param y Second line.
 * @param tol Tolerance.
 * @return true when parallel.
 */
bool c_is_parallel_lines(const c_line *x, const c_line *y, double tol);
/**
 * Test whether two planes are parallel.
 * @param x First plane.
 * @param y Second plane.
 * @param tol Tolerance.
 * @return true when parallel.
 */
bool c_is_parallel_planes(const c_plane *x, const c_plane *y, double tol);
/**
 * Test whether a point lies on a plane.
 * @param pt Point.
 * @param pln Plane.
 * @param tol Tolerance.
 * @return true when on the plane.
 */
bool c_is_point_on_plane(const double pt[3], const c_plane *pln, double tol);
/**
 * Test whether a point lies on a line.
 * @param pt Point.
 * @param ln Line.
 * @param tol Tolerance.
 * @return true when on the line.
 */
bool c_is_point_on_line(const double pt[3], const c_line *ln, double tol);
/**
 * Find the parameter of the nearest point on a line.
 * @param pt Point.
 * @param ln Line.
 * @return Line parameter.
 */
double c_nearest_point_on_line(const double pt[3], const c_line *ln);
/**
 * Compute the distance from a point to a line.
 * @param pt Point.
 * @param ln Line.
 * @return Distance.
 */
double c_point_to_line_distance(const double pt[3], const c_line *ln);
/**
 * Compute the distance from a point to a plane.
 * @param pt Point.
 * @param pln Plane.
 * @return Distance.
 */
double c_point_to_plane_distance(const double pt[3], const c_plane *pln);
/**
 * Project a vector onto a plane.
 * @param x Input vector.
 * @param pln Plane.
 * @param px Output projection.
 */
void c_vector_plane_projection(const double x[3], const c_plane *pln, double px[3]);
/**
 * Project a point onto a plane.
 * @param pt Input point.
 * @param pln Plane.
 * @param ppt Output projected point.
 */
void c_point_plane_projection(const double pt[3], const c_plane *pln, double ppt[3]);
/**
 * Construct a Pluecker line through two points.
 * @param pt1 First point.
 * @param pt2 Second point.
 * @param ln Output Pluecker line.
 */
void c_plucker_line_from_2_points(const double pt1[3], const double pt2[3], 
    c_plucker_line *ln);
/**
 * Convert a line to Pluecker coordinates.
 * @param src Input line.
 * @param ln Output Pluecker line.
 */
void c_plucker_line_from_line(const c_line *src, c_plucker_line *ln);
/**
 * Construct a Pluecker line from two planes.
 * @param p1 First plane.
 * @param p2 Second plane.
 * @param ln Output line.
 */
void c_plucker_line_from_2_planes(const c_plane *p1, const c_plane *p2, 
    c_plucker_line *ln);
/**
 * Construct a Pluecker line from six coordinates.
 * @param x Six coordinates.
 * @param ln Output line.
 */
void c_plucker_line_from_array(const double x[6], c_plucker_line *ln);
/**
 * Apply a matrix to a Pluecker line.
 * @param n Matrix row count.
 * @param x Matrix.
 * @param ldx Leading dimension of x.
 * @param ln Pluecker line.
 * @param y Output vector.
 */
void c_plucker_line_mtx_mult(int n, const double *x, int ldx, 
    const c_plucker_line *ln, double *y);
/**
 * Convert a Pluecker line to six coordinates.
 * @param ln Input line.
 * @param x Output coordinates.
 */
void c_plucker_line_to_array(const c_plucker_line *ln, double x[6]);
/**
 * Compute the common normal of two lines.
 * @param ln1 First line.
 * @param ln2 Second line.
 * @param ln Output common normal.
 */
void c_line_common_normal(const c_line *ln1, const c_line *ln2, c_line *ln);
/**
 * Determine whether two lines intersect and return their parameters.
 * @param ln1 First line.
 * @param ln2 Second line.
 * @param intersect Output intersection flag.
 * @param t1 Output parameter on ln1.
 * @param t2 Output parameter on ln2.
 * @param tol Tolerance.
 */
void c_do_lines_intersect(const c_line *ln1, const c_line *ln2, bool *intersect,
    double *t1, double *t2, double tol);
/**
 * Construct a line from a point and direction.
 * @param pt Point.
 * @param v Direction vector.
 * @param ln Output line.
 */
void c_line_from_point_and_vector(const double pt[3], const double v[3],
    c_line *ln);

/**
 * Compute a Poincare section map from sampled trajectories.
 * @param n Sample count.
 * @param x X samples.
 * @param y Y samples.
 * @param z Z samples.
 * @param pln Section plane.
 * @param side Section orientation selector.
 * @param nbuffer Buffer capacity.
 * @param xbuffer Output section x values.
 * @param ybuffer Output section y values.
 * @param zbuffer Output section z values.
 * @param nactual Output number of intersections.
 */
void c_poincare_map(int n, const double *x, const double *y, const double *z,
    const c_plane *pln, int side, int nbuffer, double *xbuffer, double *ybuffer,
    double *zbuffer, int *nactual);
/**
 * @}
 */

/**
 * @defgroup dynamics_serial Serial linkage operations
 */
/**
 * @{
 */
/**
 * Allocate a Denavit-Hartenberg table.
 * @param n Number of parameter sets.
 * @param tbl Table initialized by the routine.
 * @return Zero on success; nonzero on allocation failure.
 */
int c_alloc_dh_table(int n, c_dh_table *tbl);
/**
 * Release a Denavit-Hartenberg table.
 * @param tbl Table to release.
 */
void c_free_dh_table(c_dh_table *tbl);
/**
 * Define a link coordinate system from adjacent geometry.
 * @param xim1 Previous x direction.
 * @param zim1 Previous z direction.
 * @param zi Current z direction.
 * @param rim1 Previous origin.
 * @param ri Current origin.
 * @param csys Output coordinate system.
 */
void c_define_link_csys(const double xim1[3], const double zim1[3], 
    const double zi[3], const double rim1[3], const double ri[3],
    c_coordinate_system *csys);
/**
 * Define a coordinate system from basis vectors and origin.
 * @param i First basis vector.
 * @param j Second basis vector.
 * @param k Third basis vector.
 * @param o Origin.
 * @param csys Output coordinate system.
 */
void c_define_csys(const double i[3], const double j[3], const double k[3], 
    const double o[3], c_coordinate_system *csys);
/**
 * Build Denavit-Hartenberg parameters from coordinate systems.
 * @param n Number of systems.
 * @param csys Coordinate systems.
 * @param tbl Output table.
 */
void c_build_dh_table(int n, const c_coordinate_system *csys, c_dh_table *tbl);

/**
 * Allocate a serial linkage.
 * @param n Link count.
 * @param lnk Linkage initialized by the routine.
 * @return Zero on success; nonzero on allocation failure.
 */
int c_alloc_serial_linkage(int n, c_serial_linkage *lnk);
/**
 * Release a serial linkage.
 * @param lnk Linkage to release.
 */
void c_free_serial_linkage(c_serial_linkage *lnk);
/**
 * Build a serial linkage from link definitions.
 * @param n Link count.
 * @param links Link definitions.
 * @param linkage Output linkage.
 */
void c_build_serial_linkage(int n, const c_binary_link *links, 
    c_serial_linkage *linkage);
/**
 * Evaluate serial-link forward kinematics.
 * @param n Link count.
 * @param lnk Linkage.
 * @param q Joint variables.
 * @param T Output transform.
 * @param ldt Leading dimension of T.
 */
void c_serial_linkage_forward_kinematics(int n, const c_serial_linkage *lnk,
    const double *q, double *T, int ldt);
/**
 * Evaluate the serial-link Jacobian.
 * @param n Link count.
 * @param lnk Linkage.
 * @param q Joint variables.
 * @param jac Output Jacobian.
 * @param ldj Leading dimension of jac.
 */
void c_serial_linkage_jacobian(int n, const c_serial_linkage *lnk, 
    const double *q, double *jac, int ldj);
/**
 * Solve serial-link inverse kinematics.
 * @param n Link count.
 * @param lnk Linkage.
 * @param qo Initial joint variables.
 * @param trg Target transform.
 * @param ldt Leading dimension of trg.
 * @param q Output joint variables.
 * @param ib Output iteration statistics.
 */
void c_serial_linkage_inverse_kinematics(int n, const c_serial_linkage *lnk,
    const double *qo, const double *trg, int ldt, double *q, 
    c_iteration_behavior *ib);
/**
 * @}
 */

/**
 * @defgroup dynamics_parallel Parallel and planar linkage operations
 */
/**
 * @{
 */
/**
 * Allocate a multi-frame mechanism link.
 * @param nframes Number of link frames.
 * @param lnk Link initialized by the routine.
 * @return Zero on success; nonzero on allocation failure.
 */
int c_alloc_mechanism_link(int nframes, c_mechanism_link *lnk);
/**
 * Release a multi-frame mechanism link.
 * @param lnk Link to release.
 */
void c_free_mechanism_link(c_mechanism_link *lnk);
/**
 * Create a spatial parallel linkage mechanism.
 * @param nlinks Link count.
 * @param links Link definitions.
 * @param njoints Joint count.
 * @param joints Joint definitions.
 * @param base One-based base-link index.
 * @param effector One-based end-effector link index.
 * @param tool Tool transform.
 * @param ldt Leading dimension of tool.
 * @return Opaque mechanism handle, or NULL on failure.
 */
c_mechanism c_create_parallel_linkage(int nlinks, const c_mechanism_link *links,
    int njoints, const c_joint *joints, int base, int effector,
    const double *tool, int ldt);
/**
 * Create a planar parallel linkage mechanism.
 * @param nlinks Link count.
 * @param links Link definitions.
 * @param njoints Joint count.
 * @param joints Joint definitions.
 * @param base One-based base-link index.
 * @param effector One-based end-effector link index.
 * @param tool Tool transform.
 * @param ldt Leading dimension of tool.
 * @return Opaque mechanism handle, or NULL on failure.
 */
c_mechanism c_create_planar_linkage(int nlinks, const c_mechanism_link *links,
    int njoints, const c_joint *joints, int base, int effector,
    const double *tool, int ldt);
/**
 * Release a mechanism handle.
 * @param obj Mechanism handle.
 */
void c_free_mechanism(c_mechanism obj);
/**
 * Return the number of links in a mechanism.
 * @param obj Mechanism handle.
 * @return Link count.
 */
int c_mechanism_link_count(c_mechanism obj);
/**
 * Return the number of joints in a mechanism.
 * @param obj Mechanism handle.
 * @return Joint count.
 */
int c_mechanism_joint_count(c_mechanism obj);
/**
 * Return the number of mechanism variables.
 * @param obj Mechanism handle.
 * @return Variable count.
 */
int c_mechanism_variable_count(c_mechanism obj);
/**
 * Return the number of independent loops.
 * @param obj Mechanism handle.
 * @return Loop count.
 */
int c_mechanism_loop_count(c_mechanism obj);
/**
 * Return the number of constraint equations.
 * @param obj Mechanism handle.
 * @return Constraint count.
 */
int c_mechanism_constraint_count(c_mechanism obj);
/**
 * Return mechanism degrees of freedom.
 * @param obj Mechanism handle.
 * @return Degrees of freedom.
 */
int c_mechanism_degrees_of_freedom(c_mechanism obj);
/**
 * Return the number of actuated variables.
 * @param obj Mechanism handle.
 * @return Actuated-variable count.
 */
int c_mechanism_actuated_variable_count(c_mechanism obj);
/**
 * Return the mechanism spatial dimension.
 * @param obj Mechanism handle.
 * @return Space dimension.
 */
int c_mechanism_space_dimension(c_mechanism obj);
/**
 * Return the frame count for a link.
 * @param obj Mechanism handle.
 * @param i One-based link index.
 * @return Frame count.
 */
int c_mechanism_link_frame_count(c_mechanism obj, int i);
/**
 * Return a link-frame transform.
 * @param obj Mechanism handle.
 * @param i One-based link index.
 * @param k One-based frame index.
 * @param T Output transform.
 * @param ldt Leading dimension of T.
 */
void c_mechanism_link_frame(c_mechanism obj, int i, int k, double *T, int ldt);
/**
 * Get the current mechanism configuration.
 * @param obj Mechanism handle.
 * @param n Configuration length.
 * @param q Output configuration.
 */
void c_mechanism_get_configuration(c_mechanism obj, int n, double *q);
/**
 * Set the mechanism configuration.
 * @param obj Mechanism handle.
 * @param n Configuration length.
 * @param q Configuration values.
 */
void c_mechanism_set_configuration(c_mechanism obj, int n, const double *q);
/**
 * Evaluate a body transform at a configuration.
 * @param obj Mechanism handle.
 * @param i One-based body index.
 * @param n Configuration length.
 * @param q Configuration values.
 * @param T Output transform.
 * @param ldt Leading dimension of T.
 */
void c_mechanism_body_transform(c_mechanism obj, int i, int n, const double *q,
    double *T, int ldt);
/**
 * Evaluate the end-effector transform.
 * @param obj Mechanism handle.
 * @param n Configuration length.
 * @param q Configuration values.
 * @param T Output transform.
 * @param ldt Leading dimension of T.
 */
void c_mechanism_end_effector_transform(c_mechanism obj, int n, 
    const double *q, double *T, int ldt);
/**
 * Evaluate mechanism constraint equations.
 * @param obj Mechanism handle.
 * @param n Configuration length.
 * @param q Configuration values.
 * @param nc Constraint count.
 * @param f Output constraints.
 */
void c_mechanism_constraints(c_mechanism obj, int n, const double *q, int nc,
    double *f);
/**
 * Evaluate the mechanism constraint Jacobian.
 * @param obj Mechanism handle.
 * @param n Configuration length.
 * @param q Configuration values.
 * @param jac Output Jacobian.
 * @param ldj Leading dimension of jac.
 */
void c_mechanism_constraint_jacobian(c_mechanism obj, int n, const double *q,
    double *jac, int ldj);
/**
 * Solve for a complete mechanism configuration.
 * @param obj Mechanism handle.
 * @param na Number of actuated variables.
 * @param qa Actuated variables.
 * @param n Configuration length.
 * @param q Configuration updated by the solver.
 * @param ib Output iteration statistics.
 */
void c_mechanism_solve_configuration(c_mechanism obj, int na, const double *qa,
    int n, double *q, c_iteration_behavior *ib);
/**
 * Evaluate forward kinematics for a mechanism.
 * @param obj Mechanism handle.
 * @param na Number of actuated variables.
 * @param qa Actuated variables.
 * @param T Output end-effector transform.
 * @param ldt Leading dimension of T.
 * @param ib Output iteration statistics.
 */
void c_mechanism_forward_kinematics(c_mechanism obj, int na, const double *qa,
    double *T, int ldt, c_iteration_behavior *ib);
/**
 * Evaluate a mechanism Jacobian.
 * @param obj Mechanism handle.
 * @param na Number of actuated variables.
 * @param qa Actuated variables.
 * @param jac Output Jacobian.
 * @param ldj Leading dimension of jac.
 */
void c_mechanism_jacobian(c_mechanism obj, int na, const double *qa, 
    double *jac, int ldj);
/**
 * Solve mechanism inverse kinematics.
 * @param obj Mechanism handle.
 * @param trg Target transform.
 * @param ldt Leading dimension of trg.
 * @param na Number of actuated variables.
 * @param qa Actuated variables, updated by the solver.
 * @param ib Output iteration statistics.
 */
void c_mechanism_inverse_kinematics(c_mechanism obj, const double *trg, int ldt,
    int na, double *qa, c_iteration_behavior *ib);
/**
 * @}
 */

/**
 * @defgroup dynamics_state Transfer functions and state-space models
 */
/**
 * @{
 */
/**
 * Allocate a polynomial.
 * @param order Polynomial order.
 * @param p Polynomial initialized by the routine.
 * @return Zero on success; nonzero on allocation failure.
 */
int c_alloc_polynomial(int order, c_polynomial *p);
/**
 * Release a polynomial.
 * @param p Polynomial to release.
 */
void c_free_polynomial(c_polynomial *p);
/**
 * Allocate a transfer function.
 * @param numer_order Numerator order.
 * @param denom_order Denominator order.
 * @param tf Transfer function initialized by the routine.
 * @return Zero on success; nonzero on allocation failure.
 */
int c_alloc_transfer_function(int numer_order, int denom_order, 
    c_transfer_function *tf);
/**
 * Release a transfer function.
 * @param tf Transfer function to release.
 */
void c_free_transfer_function(c_transfer_function *tf);
/**
 * Allocate a continuous state-space model.
 * @param dimension State dimension.
 * @param n_inputs Input count.
 * @param n_outputs Output count.
 * @param mdl Model initialized by the routine.
 * @return Zero on success; nonzero on allocation failure.
 */
int c_alloc_state_space_model(int dimension, int n_inputs, int n_outputs, 
    c_state_space_model *mdl);
/**
 * Release a state-space model.
 * @param mdl Model to release.
 */
void c_free_state_space_model(c_state_space_model *mdl);
/**
 * Evaluate a transfer function at complex points.
 * @param tf Transfer function.
 * @param n Evaluation count.
 * @param s Complex evaluation points.
 * @param z Output complex values.
 */
void c_evaluate_transfer_function(const c_transfer_function *tf, int n,
    const double complex *s, double complex *z);
/**
 * Compute transfer-function poles.
 * @param tf Transfer function.
 * @param n Pole count.
 * @param p Output poles.
 */
void c_transfer_function_poles(const c_transfer_function *tf, int n, 
    double complex *p);
/**
 * Compute transfer-function zeros.
 * @param tf Transfer function.
 * @param n Zero count.
 * @param z Output zeros.
 */
void c_transfer_function_zeros(const c_transfer_function *tf, int n,
    double complex *z);
/**
 * Convert a transfer function to controllable canonical state space.
 * @param tf Transfer function.
 * @param ss Output state-space model.
 */
void c_to_ccf_state_space(const c_transfer_function *tf, c_state_space_model *ss);
/**
 * Convert a transfer function to observable canonical state space.
 * @param tf Transfer function.
 * @param ss Output state-space model.
 */
void c_to_ocf_state_space(const c_transfer_function *tf, c_state_space_model *ss);
/**
 * Create a state-space model from mass, damping, and stiffness matrices.
 * @param n State dimension.
 * @param n_out Output count.
 * @param m Mass matrix.
 * @param ldm Leading dimension of m.
 * @param b Damping matrix.
 * @param ldb Leading dimension of b.
 * @param k Stiffness matrix.
 * @param ldk Leading dimension of k.
 * @param mdl Output model.
 */
void c_create_state_space_model(int n, int n_out, const double *m, int ldm,
    const double *b, int ldb, const double *k, int ldk, 
    c_state_space_model *mdl);
/**
 * Create a PID-controlled plant state-space model.
 * @param kp Proportional gain.
 * @param ki Integral gain.
 * @param kd Derivative gain.
 * @param tau Derivative filter time constant.
 * @param plant Plant model.
 * @param mdl Output model.
 */
void c_create_pid_state_space_model(double kp, double ki, double kd, double tau,
    const c_state_space_model *plant, c_state_space_model *mdl);
/**
 * Multiply two transfer functions.
 * @param tf1 First transfer function.
 * @param tf2 Second transfer function.
 * @param tf Output product.
 */
void c_transfer_function_multiply(const c_transfer_function *tf1,
    const c_transfer_function *tf2, c_transfer_function *tf);
/**
 * Scale a transfer function.
 * @param x Scale factor.
 * @param tf1 Input transfer function.
 * @param tf Output scaled transfer function.
 */
void c_scale_transfer_function(double x, const c_transfer_function *tf1,
    c_transfer_function *tf);
/**
 * Integrate a continuous state-space model.
 * @param mdl State-space model.
 * @param u Input callback.
 * @param n Time-sample count.
 * @param t Time samples.
 * @param ndof Input degrees of freedom.
 * @param ic Initial state.
 * @param solver Integration method.
 * @param nout Output count.
 * @param y Output samples.
 * @param ldy Leading dimension of y.
 */
void c_lti_solve(const c_state_space_model *mdl, const c_ss_excitation u,
    int n, const double *t, int ndof, const double *ic, int solver, 
    int nout, double *y, int ldy);
/**
 * Compute state-space poles.
 * @param mdl State-space model.
 * @param n Pole count.
 * @param p Output poles.
 */
void c_state_space_poles(const c_state_space_model *mdl, int n,
    double complex *p);
/**
 * Compute state-space zeros.
 * @param mdl State-space model.
 * @param n Zero count.
 * @param z Output zeros.
 * @param nz Output number of zeros.
 */
void c_state_space_zeros(const c_state_space_model *mdl, int n,
    double complex *z, int *nz);
/**
 * Evaluate the state-space transfer matrix.
 * @param mdl State-space model.
 * @param nin Input count.
 * @param nout Output count.
 * @param n Evaluation count.
 * @param s Complex evaluation points.
 * @param z Output transfer values.
 * @param ldz Leading dimension of z.
 */
void c_state_space_transfer_function(const c_state_space_model *mdl, int nin,
    int nout, int n, const double complex *s, double complex *z, int ldz);
/**
 * @}
 */

#ifdef __cplusplus
}
#endif
#endif



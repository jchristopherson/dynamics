# dynamics
A library of routines used for the analysis of dynamic systems.

## Status
[![CMake](https://github.com/jchristopherson/dynamics/actions/workflows/cmake.yml/badge.svg)](https://github.com/jchristopherson/dynamics/actions/workflows/cmake.yml)
[![Actions Status](https://github.com/jchristopherson/dynamics/workflows/fpm/badge.svg)](https://github.com/jchristopherson/dynamics/actions)

## Documentation
- [Fortran API documentation](https://jchristopherson.github.io/dynamics/)
- [C API documentation](https://jchristopherson.github.io/dynamics/c-api/)

## Capabilities
The `dynamics` module aggregates tools for analysis, modeling, and identification of dynamic systems.

- Frequency response and modal analysis
    - SISO and MIMO FRF computation for linear systems.
    - Modal response and proportional/modal damping support.
    - Nonlinear frequency sweep workflows (ascending/descending) to expose effects such as jump behavior.
    - FRF model fitting (for example, accelerance and receptance models).
- Vibrations and response characterization
    - Q-factor and bandwidth estimation.
    - Damping estimation (for example, logarithmic decrement and overshoot-based methods).
    - Free-response property extraction (resonant frequency, damping ratio, settling amplitude, etc.).
    - Step-response metrics such as rise time and settling behavior.
- Controls and system representations
    - State-space and transfer-function representations.
    - LTI simulation utilities and polynomial helpers.
- System identification
    - Least-squares parameter estimation of dynamic models from measured input/output data.
    - Regression statistics and solver controls for fit quality and convergence behavior.
- Kinematics, rigid-body motion, and robotics utilities
    - Denavit-Hartenberg parameter sets, transformations, and forward/inverse kinematics. See the [DH parameter diagram](images/dh_parameter_set.svg).
    - Jacobian-related helpers for mechanism analysis.
    - Serial-link linkage modeling (including revolute/prismatic joint handling).
    - Closed-loop (parallel) mechanism modeling with loop-closure constraints, mobility calculations, and constraint-partitioned Jacobians.
    - Graph-based mechanism topology utilities (spanning trees, independent loop identification).
    - Rotation transforms, angle-axis conversion, and quaternion algebra.
- Variational multibody integration
    - Structure-preserving rigid-body integration in maximal coordinates using the formulation of Brüdigam et al. (2023).
    - Direct dynamic analysis of serial, spatial parallel, and planar parallel linkages using link mass properties and joint attachment frames.
    - World-frame joint reaction forces and moments recovered from dynamic-analysis constraint multipliers.
    - Tension/compression linear springs with free-length preload and axial-only linear viscous dampers between body or ground attachment points.
    - Revolute-joint torsional springs with free-angle preload and twist-rate-only torsional dampers, with extensible force-law base types for future nonlinear elements.
    - Holonomic equality constraints enforced at the position level with Lagrange multipliers.
    - Unit-quaternion orientation updates with body-frame angular velocities and inertia tensors.
    - Dense LU and graph-factorized block solvers for the coupled Newton equations.
    - Callback interfaces for applied forces, body-frame torques, constraints, and optional analytic constraint Jacobians.
- Geometry and vector utilities
    - Point, plane, line, and Plucker-line representations and constructors. See the [geometry operations diagram](images/geometry_operations.svg).
    - Point/line/plane projection and distance calculations.
    - Intersection/parallelism checks and common-normal calculations.
    - Vector helper routines such as cross products and skew-symmetric forms.
- Structural dynamics
    - 2D/3D beam element utilities and material/node/element abstractions, with local-coordinate [beam system documentation](images/beam_coordinate_system.svg).
    - 2D/3D axial truss elements with consistent mass and no rotational degrees of freedom.
    - Position-dependent beam shear-force and bending-moment extraction in 2D and 3D.
    - Connectivity matrix construction and boundary-condition application.
    - Sparse/CSR-oriented structural assembly helpers.
    - Generalized-alpha time stepping for linear structural systems with dense or CSR mass, damping, and stiffness matrices.
- Stability analysis
    - Local fixed-point stability classification helpers.

## Building

### Prerequisites
- A Fortran compiler with Fortran 2018 support.
- CMake 3.24+ for the CMake workflow.
- fpm for the FPM workflow.
- BLAS and LAPACK available to the linker (required by the FPM build configuration).

### Build With CMake
Configure and build from the repository root:

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

Useful CMake options:
- `-DBUILD_TESTING=ON` to build test targets.
- `-DBUILD_DYNAMICS_EXAMPLES=ON` to build example programs.
- `-DBUILD_DYNAMICS_C_INTERFACE=ON` to include the C interface.
- `-DBUILD_SHARED_LIBS=ON` to build shared libraries.

Run tests (if enabled):

```sh
ctest --test-dir build --output-on-failure
```

Install the library:

```sh
cmake --install build --prefix <install-prefix>
```

### Build With FPM
From the repository root:

```sh
fpm build --profile release
```

Run tests:

```sh
fpm test
```

Install:

```sh
fpm install --prefix <install-prefix>
```

FPM resolves the package dependencies declared in `fpm.toml` automatically.

## Quick Start

### Use From FPM
Add `dynamics` to your `fpm.toml` dependencies:

```toml
[dependencies]
dynamics = { git = "https://github.com/jchristopherson/dynamics.git", tag = "v1.4.2" }
```

Then build and run your project:

```sh
fpm build
fpm run
```

### Use From CMake
If `dynamics` is installed and discoverable via `CMAKE_PREFIX_PATH`, link it as a package:

```cmake
find_package(dynamics REQUIRED)
target_link_libraries(your_target PRIVATE dynamics::dynamics)
```

If you prefer vendoring it directly, add it as a subdirectory and link the target:

```cmake
add_subdirectory(path/to/dynamics)
target_link_libraries(your_target PRIVATE dynamics)
```

## Kinematics Example
The [`kinematics_example_1`](examples/kinematics_example_1.f90) example illustrates the forward and inverse kinematic models of the illustrated 3R mechanism. This example is Example 127 from Jazar's text "Theory of Applied Robotics, Kinematics, Dynamics, & Control."

![](images/3R%20Manipulator.PNG?raw=true)

```fortran
real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
type(binary_link), dimension(3) :: links
type(serial_linkage) :: linkage
real(real64), dimension(3) :: theta, q
real(real64), dimension(4,4) :: target

! Define the links using Denavit-Hartenberg parameters.
links(1) = binary_link(twist = 0.5d0 * pi, jtype = REVOLUTE_JOINT)
links(2) = binary_link(length = 2.0d1, offset = 1.5d0, &
    jtype = REVOLUTE_JOINT)
links(3) = binary_link(length = 1.0d1, jtype = REVOLUTE_JOINT)
linkage = serial_linkage(links)

! Compute a pose, then recover its joint variables from a zero initial guess.
call random_number(theta)
target = linkage%forward_kinematics(theta)
q = linkage%inverse_kinematics([0.0d0, 0.0d0, 0.0d0], target)
```
The output of the forward kinematics is the 4-by-4 transformation matrix relating the end-effector coordinate frame to the base coordinate frame.

```math
T = \begin{bmatrix} 0.63810402550972301 & -0.68876728703171164 & 0.34412625145905751 & 21.729664033078443 \\ 0.23387248611350506 & -0.25244115588063465 & -0.93892338508354212 & 6.3665980889762501 \\ 0.73357134136181679 & 0.67961245363267497 & 6.1230317691118863E-017 & 19.601355890978652 \\ 0 & 0 & 0 & 1 \end{bmatrix}
```

This forward result was arrived at for the following values of each joint variable (units = radians).

```math
\theta = \begin{Bmatrix} 0.35130804065430021 \\ 0.66020922397550363 \\ 0.16335289838907696 \end{Bmatrix}
```

The inverse model computed these joint variables, starting from a zero condition, as follows.

```math
\theta_{inv} = \begin{Bmatrix} 0.35130804065430021 \\ 0.66020922397550375 \\ 0.16335289838907677 \end{Bmatrix}
```

## Parallel Linkage Example
A closed-loop, or parallel, mechanism is described by the `parallel_linkage` type for spatial mechanisms and by the `planar_linkage` type for mechanisms restricted to planar motion.  The topology is supplied as a collection of `link` objects connected by `joint` objects; internally the mechanism is stored as a graph whose vertices are the links and whose edges are the joints.  A spanning tree of that graph provides the transformation path to each link, and the edges excluded from the tree define the loop-closure constraints.

Unlike a serial linkage, the forward kinematics of a closed-loop mechanism require the solution of these constraints.  As a mechanism admits more than one assembly mode, the starting estimate supplied by `set_configuration` selects the branch of interest.

The [`four_bar_example_1`](examples/four_bar_example_1.f90) example analyzes a planar four-bar linkage driven at the crank.

```fortran
type(link_container), dimension(4) :: links
type(joint), dimension(4) :: joints
type(planar_linkage) :: linkage
real(real64), dimension(100,3) :: path
real(real64), allocatable, dimension(:,:) :: jacobian

! Each link carries a joint frame at either end.
allocate(links(1)%item, source = planar_link(4.0d0)) ! ground
allocate(links(2)%item, source = planar_link(1.0d0)) ! crank
allocate(links(3)%item, source = planar_link(3.5d0)) ! coupler
allocate(links(4)%item, source = planar_link(3.0d0)) ! rocker

joints(1) = joint(REVOLUTE_JOINT, 1, 2, 1, 1, actuated = .true.)
joints(2) = joint(REVOLUTE_JOINT, 2, 3, 2, 1)
joints(3) = joint(REVOLUTE_JOINT, 3, 4, 2, 2)
joints(4) = joint(REVOLUTE_JOINT, 4, 1, 1, 2)
linkage = planar_linkage(links, joints, base = 1, effector = 3)

! Select an assembly mode, sweep the crank, and query its Jacobian.
call linkage%set_configuration([0.0d0, 0.5d0*pi, -0.5d0*pi, 0.0d0])
do i = 1, size(path,1)
    path(i,:) = linkage%end_effector_pose([2.0d0*pi*(i-1)/99.0d0])
end do
jacobian = linkage%jacobian([0.5d0*pi])
```

The program produces the following output.

```text
Number of independent loops: 1
Number of joint variables: 4
Number of constraint equations: 3
Degrees of freedom: 1

Jacobian at a crank angle of 90 degrees:
 -0.50785340914727206
 -0.10495163408386032
  0.21614053493250662
```

The mobility of the mechanism follows from the number of joint variables and the number of loop-closure constraints.

```math
\text{dof} = 4 - 3 = 1
```

The Jacobian matrix is formed by partitioning the constraint Jacobian into its actuated and passive terms.  The passive joint velocities follow from

```math
\dot{q}_{p} = -C_{p}^{-1} C_{a} \dot{q}_{a},
```

which is then combined with the end-effector Jacobian.  The partition becomes singular when the mechanism reaches a configuration in which it loses control of one or more of its degrees of freedom.

Notice that the linkage is drawn by querying the mechanism itself.  The `body_transform` routine locates the body frame of each link, and the joint frames carried by each link, available via `get_link`, then locate every joint.  This approach requires no knowledge of the geometry beyond what was used to construct the mechanism, and it extends without modification to links carrying more than two joints.

![](images/four_bar_example_1.png?raw=true)

## Variational Integrator Example
The [`variational_integrator_example`](examples/variational_integrator_example.f90) simulates a planar double pendulum in maximal coordinates. Both connecting rods have distributed mass, finite cross-section inertia, and gravity loading at their centers of mass. Six holonomic constraints pin the first rod to ground and join the two rod endpoints.

The example selects the graph-factorized solver from Brüdigam et al. (2023), supplies force and constraint callbacks, and provides an analytic reduced constraint Jacobian for efficient Newton iterations:

```fortran
type(rigid_body), dimension(2) :: bodies
type(variational_state) :: initial_state
type(variational_state), allocatable, dimension(:) :: solution
type(variational_integrator) :: integrator

! Each connecting rod carries its own mass and center-of-mass inertia tensor.
bodies(1) = rigid_body(mass1, rod_inertia(mass1, length1, width1))
bodies(2) = rigid_body(mass2, rod_inertia(mass2, length2, width2))

call initialize_variational_state(initial_state, 2)
initial_state%orientation(1) = quaternion(angle1_initial, &
    [0.0d0, 0.0d0, 1.0d0])
initial_state%orientation(2) = quaternion(angle2_initial, &
    [0.0d0, 0.0d0, 1.0d0])

! Set compatible center-of-mass positions for the two endpoint constraints.
direction1 = [sin(angle1_initial), -cos(angle1_initial), 0.0d0]
direction2 = [sin(angle2_initial), -cos(angle2_initial), 0.0d0]
initial_state%position(:,1) = 0.5d0 * length1 * direction1
initial_state%position(:,2) = length1 * direction1 + &
    0.5d0 * length2 * direction2

integrator%settings%linear_solver = VI_GRAPH_FACTORIZED_SOLVER
solution = integrator%solve(bodies, initial_state, dt, ntime, &
    constraint_count = 6, &
    constraint = pendulum_constraints, &
    force_function = gravity_forces, &
    constraint_jacobian = pendulum_constraint_jacobian, &
    args = parameters)
```
![Double-pendulum rod angles produced by the variational integrator example](images/variational_integrator_example.png?raw=true)

## Motor-Driven Parallel Linkage Example
The linkage dynamics analysis capabilities also allow the addition of motors to drive motion, along with spring and damper elements.  The [`four_bar_example_2`](examples/four_bar_example_2.f90) analyzes a four-bar linkage driven at the crank by a motor, and utilizes a torsional spring and damper at the rocker-ground revolute joint.  The analysis illustrates how, given the motor motion, to extract the motor torque required to achieve the motion, along with the joint reaction forces in terms of the world coordinate frame.  The solver always places prescribed-motion constraints after all the joint constraints; therefore, the motor torque is the last row in the Lagrange multiplier output.  For this example, the motor torque is the required torque to overcome the inertia of the mechanism, and the spring and damper at the rocker-ground revolute joint.

```fortran
type(link_container) :: links(4)
type(joint) :: joints(4)
type(planar_linkage) :: mechanism
type(linkage_dynamic_model) :: model
procedure(linkage_prescribed_motion), pointer :: motor_motion
type(torsional_spring) :: spring
type(torsional_damper) :: damper
real(real64) :: g(3), R(3,3)
type(variational_integrator) :: integrator
type(variational_state), allocatable, dimension(:) :: solution
real(real64), allocatable, dimension(:,:) :: constraint_multipliers
type(joint_reaction) :: joint_reactions(4)

! Connect the links via revolute joints
joints(1) = joint( &
    REVOLUTE_JOINT, &       ! joint type
    1, &                    ! parent link index (ground link)
    2, &                    ! child link index (crank link)
    1, &                    ! parent coordinate frame index (interfaces to joint #1 on the ground link)
    1, &                    ! child coordinate frame index  (interfaces to joint #1 on the crank link)
    actuated = .true. &     ! we're going to drive this joint with a motor
)
joints(2) = joint( &
    REVOLUTE_JOINT, &       ! joint type
    2, &                    ! parent link index (crank link)
    3, &                    ! child link index (coupler link)
    2, &                    ! parent coordinate frame index (interfaces to joint #2 on the crank)
    1 &                     ! child coordinate frame index (interfaces to joint #1 on the coupler)
)
joints(3) = joint( &
    REVOLUTE_JOINT, &       ! joint type
    3, &                    ! parent link index (coupler link)
    4, &                    ! child link index (rocker link)
    2, &                    ! parent coordinate frame index (interfaces to joint #2 on the coupler)
    2 &                     ! child coordinate frame index (interfaces to joint #2 on the rocker)
)
joints(4) = joint( &
    REVOLUTE_JOINT, &       ! joint type
    4, &                    ! parent link index (rocker link)
    1, &                    ! child link index (ground link)
    1, &                    ! parent coordinate frame index (interfaces to joint #1 on the rocker)
    2 &                     ! child coordinate frame index (interfaces to joint #2 on the ground link)
)

! Construct the mechanism
mechanism = planar_linkage( &
    links, &                ! the list of link objects
    joints, &               ! the list of joint objects
    base = 1, &             ! the index of the base link (ground)
    effector = 3 &          ! the index of the end-effector link (coupler)
)

! Define the dynamic model
model = linkage_dynamic_model( &
    mechanism, &                        ! the mechanism
    mechanism%get_configuration() &     ! the initial configuration
)

! Now we can define the necessary routines for the integrator
motor_motion => crank_motor

! Add a spring and damper element
spring%free_angle = 0.0d0
spring%joint_index = 4      ! the spring is tied to the rocker-ground connection
spring%stiffness = torsional_spring_rate
call model%add_torsional_spring(spring)

damper%joint_index = 4      ! the damper is tied to the rocker-ground connection
damper%damping = torsional_damping_rate
call model%add_torsional_damper(damper)

! Define the gravitational vector
g = [0.0d0, -gc, 0.0d0]     ! gravity in -y direction

! Set up the integrator & solve
integrator%settings%linear_solver = VI_DENSE_SOLVER ! define the solver type
solution = model%solve( &
    integrator, &                           ! the integrator to utilize
    dt, &                                   ! time step size
    ntime, &                                ! # of time steps
    gravity = g, &                          ! gravitational vector
    prescribed_body = 1, &                  ! the index of the driven body - ground is ignored in this instance, so the driving body is the crank and thus index 1
    prescribed_motion = motor_motion, &     ! motor motion routine
    multipliers = constraint_multipliers &  ! Lagrange multiplier values - use to get motor torque
)

! Extract the solution components
do i = 1, ntime
    time(i) = solution(i)%time

    ! Crank Angle
    R = solution(i)%orientation(1)%to_matrix()
    crank_angle(i) = atan2(R(2,1), R(1,1)) * 1.8d2 / pi

    ! Rocker Pivot Angle
    R = solution(i)%orientation(3)%to_matrix()
    rocker_angle(i) = atan2(R(2,1), R(1,1)) * 1.8d2 / pi

    ! Joint Reactions
    joint_reactions = model%get_joint_reactions( &
        solution(i), &
        constraint_multipliers(:,i) &
    )
    
    Fx1(i) = joint_reactions(1)%force(1)
    Fy1(i) = joint_reactions(1)%force(2)

    Fx2(i) = joint_reactions(2)%force(1)
    Fy2(i) = joint_reactions(2)%force(2)

    Fx3(i) = joint_reactions(3)%force(1)
    Fy3(i) = joint_reactions(3)%force(2)

    Fx4(i) = joint_reactions(4)%force(1)
    Fy4(i) = joint_reactions(4)%force(2)
end do

! Get the motor torque
nmult = size(constraint_multipliers, 1)
torque = constraint_multipliers(nmult,:)    ! motor torque
```
![](images/four_bar_example_2a.png?raw=true)
![](images/four_bar_example_2b.png?raw=true)

### Nonconservative Loads and Damping
The variational integrator accepts applied forces and torques, including dissipative loads such as viscous dampers. They enter as nonconservative forces in the forced discrete Euler-Lagrange equations. The `force_evaluation` setting selects when **all** applied loads are sampled; it changes the force quadrature, not the conservative state update or its formal order.

Let $t_k,y_k$ be the current time and state, and $t_{k+1},y_{k+1}$ the trial next time and state. The translational momentum balance uses the selected force sample $F_*$:
$$
m\frac{v_{k+1}-v_k}{h}=F_*.
$$
The rotational discrete momentum balance uses the torque sampled at the same point. The available modes are:

| Mode | Force sample | Tradeoffs |
|---|---|---|
| `VI_FORCE_LEFT_ENDPOINT` (default) | $F(t_k,y_k)$ | Explicit in the step and least expensive. Large damping or a large step can make the discrete velocity grow. |
| `VI_FORCE_IMPLICIT_ENDPOINT` | $F(t_{k+1},y_{k+1})$ | Evaluates loads on each Newton trial. More robust for stiff damping, but adds nonlinear work and endpoint-samples every load, including springs and user-defined loads. |
| `VI_FORCE_MIDPOINT` | $F(t_{k+1/2},y_{k+1/2})$ | Evaluates loads on each Newton trial at an interpolated midpoint. A centered force sample, but not a complete implicit-midpoint integrator and not L-stable for very stiff modes. |

For `VI_FORCE_MIDPOINT`, time, position, and translational velocity are arithmetic midpoints. Orientation follows the unit-quaternion geodesic midpoint
$$
q_{k+1/2}=q_k\left(q_k^{-1}q_{k+1}\right)^{1/2}.
$$
Angular velocity is averaged in world coordinates, then expressed in the midpoint body frame:
$$
\omega_{k+1/2}=R(q_{k+1/2})^T\frac{R(q_k)\omega_k+R(q_{k+1})\omega_{k+1}}{2}.
$$

For an isolated mass $m$ with viscous force $F=-cv$, define $r=hc/m$. The velocity amplification factors are
$$
\frac{v_{k+1}}{v_k}=\begin{cases}
1-r, & \text{left endpoint},\\
\dfrac{1}{1+r}, & \text{implicit endpoint},\\
\dfrac{1-r/2}{1+r/2}, & \text{midpoint}.
\end{cases}
$$
The left-endpoint update can grow kinetic energy when $r>2$. The implicit endpoint is unconditionally stable for this scalar linear problem and strongly damps stiff modes. The midpoint factor does not grow for $r\ge0$, but approaches $-1$ as damping becomes very stiff, so those modes alternate rather than being rapidly suppressed. These scalar properties do not guarantee an energy law for a constrained multibody simulation; check timestep convergence and energy/work behavior for the system of interest.

Select a Fortran mode with `integrator%settings%force_evaluation`, for example `VI_FORCE_MIDPOINT`. In C, set `c_variational_integrator_settings.force_evaluation` to `DYN_VI_FORCE_LEFT_ENDPOINT`, `DYN_VI_FORCE_IMPLICIT_ENDPOINT`, or `DYN_VI_FORCE_MIDPOINT`; initialize the settings with `c_default_variational_integrator_settings` first. Force callbacks may be reevaluated repeatedly during Newton iterations and should compute loads deterministically from their input time, state, and user data.

The complete example includes the massive-rod inertia calculation, gravity and endpoint-constraint callbacks, analytic quaternion-tangent Jacobian, and plots of both rod angles. Build it with `BUILD_DYNAMICS_EXAMPLES=ON` and run the `variational_integrator_example` target.

## Frequency Response Example
Consider the following 3 DOF system. The [`frf_proportional_example_1`](examples/frf_proportional_example_1.f90) example illustrates how to use this library to compute the frequency response functions for this system.

![](images/3%20DOF%20Schematic.PNG?raw=true)

The equations describing this system are as follows.

```math
\begin{bmatrix} m_1 & 0 & 0 \\ 0 & m_2 & 0 \\ 0 & 0 & m_3 \end{bmatrix} \begin{Bmatrix} \ddot{x}_1 \\ \ddot{x}_2 \\ \ddot{x}_3 \end{Bmatrix} + \begin{bmatrix} b_1 + b_2 & -b_2 & 0 \\ -b_2 & b_2 + b_3 & -b_3 \\ 0 & -b_3 & b_3 + b_4 \end{bmatrix} \begin{Bmatrix} \dot{x}_1 \\ \dot{x}_2 \\ \dot{x}_3 \end{Bmatrix} + \begin{bmatrix} k_1 + k_2 & -k_2 & 0 \\ -k_2 & k_2 + k_3 & -k_3 \\ 0 & -k_3 & k_3 + k_4 \end{bmatrix} \begin{Bmatrix} x_{1} \\ x_{2} \\ x_{3} \end{Bmatrix} = \begin{Bmatrix} F(t) \\ 0 \\ 0 \end{Bmatrix}
```

This analysis makes use of proportional damping.  Using proportional damping, the damping matrix is determined as follows.

```math
B = \alpha M + \beta K
```

The essential excitation and solution setup is:

```fortran
real(real64), dimension(3,3) :: mass, stiffness
type(frf) :: response
procedure(modal_excite), pointer :: excitation

mass = reshape([0.5d0, 0.0d0, 0.0d0, &
    0.0d0, 2.5d0, 0.0d0, &
    0.0d0, 0.0d0, 0.75d0], [3,3])
stiffness = reshape([15.0d6, -10.0d6, 0.0d0, &
    -10.0d6, 20.0d6, -10.0d6, &
    0.0d0, -10.0d6, 15.0d6], [3,3])

excitation => modal_frf_forcing_term
response = frequency_response(mass, stiffness, 1.0d-3, 2.0d-6, &
    1000, 2.0d0*pi*10.0d0, 2.0d0*pi*1.0d3, excitation)

contains
subroutine modal_frf_forcing_term(freq, force, args)
    real(real64), intent(in) :: freq
    complex(real64), intent(out), dimension(:) :: force
    class(*), intent(inout), optional :: args
    force = [(1.0d3, 0.0d0), (0.0d0, 0.0d0), (0.0d0, 0.0d0)]
end subroutine
```

The computed frequency response functions.

![](images/frf_proportional_example_1.png?raw=true)

## Nonlinear FRF Example
Computing the frequency response function for a nonlinear system is not as straight-forward. A technique for capturing nonlinear behaviors, such as jump phenomenon, is to sweep through frequency, in both an ascending and a descending manner. The [`frf_sweep_example_1`](examples/frf_sweep_example_1.f90) example illustrates such a frequency sweep using the famous Duffing equation as the model.

```math
\ddot{x} + \delta \dot{x} + \alpha x + \beta x^3 = \gamma \sin \omega t
```

The essential model and sweep calls are:

```fortran
procedure(harmonic_ode), pointer :: model
type(frf) :: ascending, descending

model => duffing_ode
ascending = frequency_sweep(model, 100, 0.5d0, 2.0d0, &
    [0.0d0, 0.0d0])
descending = frequency_sweep(model, 100, 2.0d0, 0.5d0, &
    [0.0d0, 0.0d0])

contains
pure subroutine duffing_ode(frequency, t, state, derivative, args)
    real(real64), intent(in) :: frequency, t
    real(real64), intent(in), dimension(:) :: state
    real(real64), intent(out), dimension(:) :: derivative
    class(*), intent(inout), optional :: args

    derivative(1) = state(2)
    derivative(2) = sin(frequency*t) - 0.1d0*state(2) - &
        state(1) - 0.04d0*state(1)**3
end subroutine
```
The computed frequency response functions, both ascending and descending, as compared with the analytical approximation.

![](images/frf_sweep_example_1.png?raw=true)

## Parameter Discovery (System Identification):
The [`siso_lsq_fit_example`](examples/siso_lsq_fit_example.f90) example illustrates how to estimate parameters of an ODE given an observed output to a known input. It finds $\omega_{n}$ and $\zeta$ in the model of a single degree of freedom system.
```math
\ddot{x} + 2 \zeta \omega_{n} \dot{x} + \omega_{n}^{2} x = f(t)
```

```fortran
type(dynamic_system_measurement), dimension(1) :: measurements
type(regression_statistics), dimension(2) :: statistics
type(iteration_controls) :: controls
procedure(ode), pointer :: model
real(real64), dimension(2) :: initial_state, parameters

! Populate measured time, input, and output arrays.
allocate(measurements(1)%t(npts), measurements(1)%input(npts), &
    measurements(1)%output(npts))
measurements(1)%t = [(dt*i, i=0, npts-1)]
measurements(1)%input = applied_force
measurements(1)%output = measured_response

! Fit natural frequency and damping ratio from an initial estimate.
parameters = [2.5d2, 1.0d-1]
initial_state = 0.0d0
model => eom
call controls%set_to_default()
call siso_model_fit_least_squares(model, measurements, initial_state, &
    parameters, controls = controls, stats = statistics)

contains
subroutine eom(t, state, derivative, args)
    real(real64), intent(in) :: t
    real(real64), intent(in), dimension(:) :: state
    real(real64), intent(out), dimension(:) :: derivative
    class(*), intent(inout), optional :: args
    ! Extract wn, zeta, and forcing from args; then evaluate the SDOF model.
end subroutine
```
The results are as follows.
```txt
NATURAL FREQUENCY TERM:
        Actual:  300.000 rad/s
        Computed:  299.974 rad/s
        Difference:   -0.026 rad/s
        Std. Error:    0.154 rad/s
        Conf. Int.: +/-   0.302 rad/s
        P-Value:  0.000E+00
        T-Statistic:    1.952E+03
DAMPING TERM:
        Actual:  0.050
        Computed:  0.049
        Difference: -0.001
        Std. Error:  0.001
        Conf. Int.: +/- 0.001
        P-Value:  0.000E+00
        T-Statistic:   68.707E+00
```
![](images/siso_least_squares_fit_example.png?raw=true)

## Structural Example

TO DO: Place example here

## Harmonic Truss Example

The [`harmonic_truss_example`](examples/harmonic_truss_example.f90) models a four-node, five-bar pin-jointed truss with a pinned support and a roller. It assembles axial stiffness and consistent translational mass matrices, applies the support constraints, and advances the response to a vertical 20 Hz load at the apex with `dense_generalized_alpha_integrator`. Fplot saves vertical apex and midspan motion and horizontal midspan and roller motion to `harmonic_truss_response.png` in the process's working directory. Build the `harmonic_truss_example` CMake target with examples enabled, then run it from the build's examples directory.

After assembling and reducing the mass and stiffness matrices, the time-stepping core is:

```fortran
call integrator%initialize(reduced_mass, damping, reduced_stiffness, rho_infinity = 0.7d0)
do time_index = 1, nsteps + 1
    time(time_index) = (time_index - 1) * dt
    apex_motion(time_index) = 1.0d3 * displacement(apex_y)
    middle_motion(time_index) = 1.0d3 * displacement(middle_y)
    if (time_index > nsteps) exit

    force_next = 0.0d0
    force_next(apex_y) = -force_amplitude * &
        sin(2.0d0 * pi * frequency * real(time_index, real64) * dt)
    call integrator%step(force_current, force_next, dt, &
        displacement, velocity, acceleration)
    force_current = force_next
end do
```

![Forced response of the simple truss example](images/harmonic_truss_example.png?raw=true)

## References
1. J. D. Hartog, "Mechanical Vibrations," New York: Dover Publications, Inc., 1985.
2. S. S. Rau, "Mechanical Vibrations," 3rd ed., Reading, MA: Addison-Wesley Publishing Co., 1995.
3. R. N. Jazar, "Advanced Vibrations," 2nd ed., New York: Springer, 2022.
4. W. T. Thomson, "Theory of Vibration with Applications," 4th ed., New York: Springer, 1993.
5. A. H. Nayfeh and B. Balachandran, "Applied Nonlinear Dynamics. Analytical, Computational, and Experimental Methods," New York: John WIley & Sons, Inc., 1995.
6. L. Meirovitch, "Fundamentals of Vibrations," Long Grove, IL: Waveland Press, Inc., 2001.
7. R. N. Jazar, "Theory of Applied Robotics, Kinematics, Dynamics, and Control," New York: Springer, 2007.
8. A. H. Nayfeh, "Introduction to Perturbation Techniques," New York: John Wiley & Sons, Inc., 1993.
9. Jolicoeur, M.P., Roumy, J.G., Vanreusel, S., Dionne, D., Douville, H., Boulet, B., Michalska, H., Masson, P., & Berry, A. (2005). "Reduction of structure-borne noise in automobiles by multivariable feedback." 1397 - 1402. 10.1109/CCA.2005.1507327. 
10. Brunton, Steven & Proctor, Joshua & Kutz, J.. (2015). "Discovering governing equations from data: Sparse identification of nonlinear dynamical systems." Proceedings of the National Academy of Sciences. 113. 3932–3937. 10.1073/pnas.1517384113. 
11. Brüdigam, Jan & Sosnowski, Stefan & Manchester, Zac & Hirche, Sandra. (2023). Variational integrators and graph-based solvers for multibody dynamics in maximal coordinates. Multibody System Dynamics. 61. 1-34. 10.1007/s11044-023-09949-x. 

## Rigid-Body and Mechanism Modules

The rigid-body functionality is split among modules for spatial transformations, serial and closed-loop kinematics, and linkage dynamics:

* **`dynamics_rigid_bodies`** defines rigid-body mass properties: mass, center-of-gravity location in body coordinates, and inertia tensor about the center of gravity.
* **`dynamics_quaternions`** provides quaternion arithmetic and conversion between quaternions, rotation matrices, angle-axis, and roll-pitch-yaw representations. Unit quaternions represent 3D orientation.
* **`dynamics_rotation`** provides rotation and homogeneous transformation matrices, translations, and velocity/acceleration coordinate transforms.
* **`dynamics_kinematics`** provides Denavit-Hartenberg transformations and tables, forward and inverse kinematics, and kinematic Jacobians.
* **`dynamics_linkage`** defines rigid links with joint attachment frames and serial linkage kinematics.
* **`dynamics_joints`** defines fixed, revolute, prismatic, cylindrical, universal, and spherical joints, including their relative degrees of freedom and constraint counts.
* **`dynamics_graph`** supplies graph, path, loop, and spanning-tree types used to represent linkage topology.
* **`dynamics_parallel_linkage`** represents closed-loop spatial and planar mechanisms, finds loop-closure constraints from the topology, and provides configuration, pose, and Jacobian queries.
* **`dynamics_linkage_dynamics`** adapts serial and planar linkages to the variational integrator, with prescribed motion, joint reactions, and axial or torsional spring/damper elements.

## Kinematics Example
The [`kinematics_example_1`](../../examples/kinematics_example_1.f90) example illustrates the forward and inverse kinematic models of the illustrated 3R mechanism. This example is Example 127 from Jazar's text "Theory of Applied Robotics, Kinematics, Dynamics, & Control."

![](../../images/3R%20Manipulator.PNG?raw=true)

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
The output of the forward kinematics is the 4-by-4 transformation matrix relating the end-effector coordinate frame to the base coordinate frame. The joint variables are generated with `random_number`, so the displayed matrix and angles are one sample; they can differ on another run. Inverse kinematics starts from the supplied zero estimate and seeks a configuration that reaches the target pose.

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

The [`four_bar_example_1`](../../examples/four_bar_example_1.f90) example analyzes a planar four-bar linkage driven at the crank.

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

![](../../images/four_bar_example_1.png?raw=true)

## Motor-Driven Parallel Linkage Example
The linkage dynamics analysis capabilities also allow the addition of motors to drive motion, along with spring and damper elements. The [`four_bar_example_2`](../../examples/four_bar_example_2.f90) analyzes a four-bar linkage driven at the crank by a motor, and utilizes a torsional spring and damper at the rocker-ground revolute joint. The analysis illustrates how, given the motor motion, to extract the motor torque required to achieve the motion, along with the joint reaction forces in terms of the world coordinate frame. The solver always places prescribed-motion constraints after all the joint constraints; therefore, the motor torque is the last row in the Lagrange multiplier output. For this example, the motor torque is the required torque to overcome the inertia of the mechanism, and the spring and damper at the rocker-ground revolute joint.

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
![](../../images/four_bar_example_2a.png?raw=true)
![](../../images/four_bar_example_2b.png?raw=true)

### Initial Discrete Momenta
`integrator%solve` interprets a new state's velocities as physical initial velocities. Its first step uses a discrete Legendre transform rather than treating those velocities as incoming discrete interval velocities. For translation, the startup balance is
$$
2m\frac{v_1-v_0}{h}=F_*+J_0^T\lambda_0.
$$
The rotational startup likewise matches the physical body-frame momentum $I\omega_0$ to the outgoing rotational discrete momentum, with half-step applied and constraint impulses. This avoids the half-acceleration startup load artifact for accelerating prescribed motion. It does not remove finite-step load errors or physical discontinuities in the prescribed motion.

The initial trajectory sample retains its physical velocities. Subsequent samples contain discrete interval velocities and have `state%discrete_momenta = .true.`. Passing a returned state into another `solve` continues without repeating startup. Reset this flag to `.false.` when replacing velocities with new physical initial values; set it to `.true.` when supplying already initialized discrete velocities. Standalone `step` calls retain their existing discrete update by default; use `initialize_momenta = .true.` on the first call to start from physical velocities. Failed startup steps leave the state and its flag unchanged.

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
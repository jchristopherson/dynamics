# dynamics
A library of routines used for the analysis of dynamic systems.

## Status
[![CMake](https://github.com/jchristopherson/dynamics/actions/workflows/cmake.yml/badge.svg)](https://github.com/jchristopherson/dynamics/actions/workflows/cmake.yml)
[![Actions Status](https://github.com/jchristopherson/dynamics/workflows/fpm/badge.svg)](https://github.com/jchristopherson/dynamics/actions)

## Documentation
The documentation can be found [here](https://jchristopherson.github.io/dynamics/).

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
    - Denavit-Hartenberg tools and forward/inverse kinematics.
    - Jacobian-related helpers for mechanism analysis.
    - Serial-link linkage modeling (including revolute/prismatic joint handling).
    - Closed-loop (parallel) mechanism modeling with loop-closure constraints, mobility calculations, and constraint-partitioned Jacobians.
    - Graph-based mechanism topology utilities (spanning trees, independent loop identification).
    - Rotation transforms, angle-axis conversion, and quaternion algebra.
- Geometry and vector utilities
    - Plane/line/plucker-line utilities.
    - Point/line/plane projection and distance calculations.
    - Intersection/parallelism checks and common-normal calculations.
    - Vector helper routines such as cross products and skew-symmetric forms.
- Structural dynamics
    - 2D/3D beam element utilities and material/node/element abstractions.
    - Connectivity matrix construction and boundary-condition application.
    - Sparse/CSR-oriented structural assembly helpers.
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
The following example illustrates the forward and inverse kinematic models of the illustrated 3R mechanism.  This example is Example 127 from Jazar's text "Theory of Applied Robotics, Kinematics, Dynamics, & Control."

![](images/3R%20Manipulator.PNG?raw=true)

```fortran
program example
    use iso_fortran_env
    use dynamics
    implicit none
    
    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Model Properties
    real(real64), parameter :: L1 = 1.5d0
    real(real64), parameter :: L2 = 2.0d1
    real(real64), parameter :: L3 = 1.0d1

    ! Local Variables
    type(binary_link) :: links(3)
    type(serial_linkage) :: linkage
    integer(int32) :: i
    real(real64) :: theta(3), T(4, 4), qo(3), q(3)

    ! Link 1 Definition
    links(1) = binary_link(twist = 0.5d0 * pi, jtype = REVOLUTE_JOINT)
    
    ! Link 2 Definition
    links(2) = binary_link(length = L2, offset = L1, jtype = REVOLUTE_JOINT)

    ! Link 3 Definition
    links(3) = binary_link(length = L3, jtype = REVOLUTE_JOINT)

    ! Build the linkage
    linkage = serial_linkage(links)

    ! Define the joint variables
    call random_number(theta)

    ! --------------------
    ! Compute the forward kinematics & display the matrix
    T = linkage%forward_kinematics(theta)
    do i = 1, 4
        print *, T(i,:)
    end do

    ! --------------------
    ! Solve the inverse problem using the end-effector position and orientation
    ! computed by the forward kinematics process as a target for the inverse
    ! calculations

    ! Define an initial guess for the joint variables
    qo = [0.0d0, 0.0d0, 0.0d0]

    ! Solve the model
    q = linkage%inverse_kinematics(qo, T)

    ! Display the solution and compare with the actual
    print "(A)", "COMPUTED (Inverse Model):"
    print *, q
    print "(A)", "ACTUAL:"
    print *, theta
end program
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

The following example analyzes a planar four-bar linkage driven at the crank.

```fortran
program example
    use iso_fortran_env
    use dynamics
    use fplot_core
    use linalg, only : identity
    implicit none

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    integer(int32), parameter :: npts = 100

    ! Model Properties
    real(real64), parameter :: crank = 1.0d0
    real(real64), parameter :: coupler = 3.5d0
    real(real64), parameter :: rocker = 3.0d0
    real(real64), parameter :: ground = 4.0d0

    ! Local Variables
    integer(int32) :: i
    type(link_container) :: links(4)
    type(joint) :: joints(4)
    type(planar_linkage) :: linkage
    real(real64) :: theta(npts), p(npts,3)
    real(real64), allocatable, dimension(:) :: plt_config
    real(real64), allocatable, dimension(:,:) :: jac

    ! Define the links.  Each link carries a joint frame at each end, with the
    ! body frame located at the first of the two.
    allocate(links(1)%item, source = planar_link(ground))
    allocate(links(2)%item, source = planar_link(crank))
    allocate(links(3)%item, source = planar_link(coupler))
    allocate(links(4)%item, source = planar_link(rocker))

    ! Connect the links.  The crank is the driven member.
    joints(1) = joint(REVOLUTE_JOINT, 1, 2, 1, 1, actuated = .true.)
    joints(2) = joint(REVOLUTE_JOINT, 2, 3, 2, 1)
    joints(3) = joint(REVOLUTE_JOINT, 3, 4, 2, 2)
    joints(4) = joint(REVOLUTE_JOINT, 4, 1, 1, 2)

    ! Build the mechanism.  The end-effector is the distal end of the coupler.
    linkage = planar_linkage(links, joints, base = 1, effector = 3, &
        tool = translation(coupler))

    ! Display the mobility of the mechanism
    print "(A,I0)", "Number of independent loops: ", linkage%get_loop_count()
    print "(A,I0)", "Number of joint variables: ", linkage%get_variable_count()
    print "(A,I0)", "Number of constraint equations: ", &
        linkage%get_constraint_count()
    print "(A,I0)", "Degrees of freedom: ", linkage%get_degrees_of_freedom()

    ! Establish a starting configuration.  A closed-loop mechanism admits more
    ! than one assembly mode, so the starting estimate selects the branch.
    call linkage%set_configuration([0.0d0, 0.5d0 * pi, -0.5d0 * pi, 0.0d0])

    ! Sweep the crank and trace the coupler point
    theta = linspace(0.0d0, 2.0d0 * pi, npts)
    do i = 1, npts
        p(i,:) = linkage%end_effector_pose([theta(i)])
    end do

    ! Examine the sensitivity of the coupler point to the crank at mid-stroke
    jac = linkage%jacobian([0.5d0 * pi])
    print "(A)", new_line('a') // "Jacobian at a crank angle of 90 degrees:"
    do i = 1, size(jac, 1)
        print *, jac(i,:)
    end do

    ! Plot the coupler curve along with the linkage. The complete set of joint 
    ! variables, not just the actuated variable, is required in order to locate 
    ! each link.
    plt_config = linkage%solve_configuration([0.25d0 * pi])
    call plot_coupler_path(p, linkage, plt_config)

contains
    ! Constructs a planar link of the requested length carrying a joint frame at
    ! each end.
    function planar_link(length) result(rst)
        real(real64), intent(in) :: length
        type(multi_joint_link) :: rst
        real(real64) :: frames(4, 4, 2)
        frames(:,:,1) = translation(0.0d0)
        frames(:,:,2) = translation(length)
        rst = multi_joint_link(frames)
    end function

    ! A 4-by-4 translation along the x-axis.
    function translation(a) result(rst)
        real(real64), intent(in) :: a
        real(real64) :: rst(4, 4)
        rst = identity(4)
        rst(1,4) = a
    end function

    ! Coupler Path Plot
    subroutine plot_coupler_path(pts, mech, q)
        real(real64), intent(in), dimension(:,:) :: pts
        class(kinematic_mechanism), intent(in), target :: mech
        real(real64), intent(in), dimension(:) :: q

        type(plot_2d) :: plt
        type(legend), pointer :: lgnd

        call plt%initialize()
        call plt%set_x_axis_title("x")
        call plt%set_y_axis_title("y")
        lgnd => plt%get_legend()
        call lgnd%set_is_visible(.true.)
        call lgnd%set_horizontal_position(LEGEND_LEFT)
        call lgnd%set_vertical_position(LEGEND_BOTTOM)
        call lgnd%set_draw_border(.false.)
        call lgnd%set_layout(LEGEND_ARRANGE_HORIZONTALLY)
        call lgnd%set_draw_inside_axes(.false.)

        ! Coupler Curve
        call plt%push(pts(:,1), pts(:,2), name = "Coupler Curve")

        ! The linkage in its requested configuration
        call draw_linkage(plt, mech, q)

        call plt%draw()
    end subroutine

    ! Draws each link as a polyline passing through the joint frames the link
    ! carries.  The location of each frame is found by transforming the frame,
    ! which is expressed in the link's body coordinate frame, into the frame of
    ! the base link.
    subroutine draw_linkage(plt, mech, q)
        type(plot_2d), intent(inout) :: plt
        class(kinematic_mechanism), intent(in), target :: mech
        real(real64), intent(in), dimension(:) :: q

        integer(int32) :: j, k, n
        class(link), pointer :: lnk
        type(plot_data_2d), allocatable, dimension(:) :: pd
        real(real64) :: T(4, 4), Pk(4, 4)
        real(real64), allocatable, dimension(:) :: x, y
        character(len = :), allocatable :: name

        allocate(pd(mech%get_link_count()))
        do j = 1, mech%get_link_count()
            ! Get the link and compute its location and orientation
            lnk => mech%get_link(j)
            n = lnk%get_joint_count()
            T = mech%body_transform(j, q)
            allocate(x(n), y(n))
            do k = 1, n
                Pk = matmul(T, lnk%get_joint_frame(k))
                x(k) = Pk(1,4)
                y(k) = Pk(2,4)
            end do

            ! Define a name for the link
            select case (j)
            case (1)
                name = "Ground Link"
            case (2)
                name = "Crank"
            case (3)
                name = "Coupler"
            case (4)
                name = "Rocker"
            end select

            ! Plot the link
            call pd(j)%define_data(x, y)
            call pd(j)%set_line_width(3.0)
            if (j == 1) call pd(j)%set_line_style(LINE_DASHED)
            call pd(j)%set_draw_markers(.true.)
            call pd(j)%set_marker_style(MARKER_EMPTY_CIRCLE)
            call pd(j)%set_name(name)
            call plt%push(pd(j))

            deallocate(x, y)
        end do
    end subroutine
end program
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

## Frequency Response Example
Consider the following 3 DOF system.  The following example illustrates how to use this library to compute the frequency response functions for this system.  

![](images/3%20DOF%20Schematic.PNG?raw=true)

The equations describing this system are as follows.

```math
\begin{bmatrix} m_1 & 0 & 0 \\ 0 & m_2 & 0 \\ 0 & 0 & m_3 \end{bmatrix} \begin{Bmatrix} \ddot{x}_1 \\ \ddot{x}_2 \\ \ddot{x}_3 \end{Bmatrix} + \begin{bmatrix} b_1 + b_2 & -b_2 & 0 \\ -b_2 & b_2 + b_3 & -b_3 \\ 0 & -b_3 & b_3 + b_4 \end{bmatrix} \begin{Bmatrix} \dot{x}_1 \\ \dot{x}_2 \\ \dot{x}_3 \end{Bmatrix} + \begin{bmatrix} k_1 + k_2 & -k_2 & 0 \\ -k_2 & k_2 + k_3 & -k_3 \\ 0 & -k_3 & k_3 + k_4 \end{bmatrix} \begin{Bmatrix} x_{1} \\ x_{2} \\ x_{3} \end{Bmatrix} = \begin{Bmatrix} F(t) \\ 0 \\ 0 \end{Bmatrix}
```

This analysis makes use of proportional damping.  Using proportional damping, the damping matrix is determined as follows.

```math
B = \alpha M + \beta K
```

The following module contains the forcing term.
```fortran
module excitation
    use iso_fortran_env
    implicit none

contains
    subroutine modal_frf_forcing_term(freq, f)
        real(real64), intent(in) :: freq
        complex(real64), intent(out), dimension(:) :: f

        complex(real64), parameter :: zero = (0.0d0, 0.0d0)
        complex(real64), parameter :: one = (1.0d0, 0.0d0)

        f = [1.0d3 * one, zero, zero]
    end subroutine
end module
```
The calling program is as follows.
```fortran
program example
    use iso_fortran_env
    use dynamics
    use excitation
    implicit none

    ! Parameters
    integer(int32), parameter :: nfreq = 1000
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: fmin = 2.0d0 * pi * 10.0d0
    real(real64), parameter :: fmax = 2.0d0 * pi * 1.0d3
    real(real64), parameter :: alpha = 1.0d-3
    real(real64), parameter :: beta = 2.0d-6

    ! Define the model parameters
    real(real64), parameter :: m1 = 0.5d0
    real(real64), parameter :: m2 = 2.5d0
    real(real64), parameter :: m3 = 0.75d0
    real(real64), parameter :: k1 = 5.0d6
    real(real64), parameter :: k2 = 10.0d6
    real(real64), parameter :: k3 = 10.0d6
    real(real64), parameter :: k4 = 5.0d6

    ! Local Variables
    real(real64) :: m(3,3), k(3,3)
    type(frf) :: rsp
    procedure(modal_excite), pointer :: fcn

    ! Initialization
    fcn => modal_frf_forcing_term

    ! Define the mass matrix
    m = reshape([m1, 0.0d0, 0.0d0, 0.0d0, m2, 0.0d0, 0.0d0, 0.0d0, m3], [3, 3])

    ! Define the stiffness matrix
    k = reshape([k1 + k2, -k2, 0.0d0, -k2, k2 + k3, -k3, 0.0d0, -k3, k3 + k4], &
        [3, 3])

    ! Compute the frequency response functions
    rsp = frequency_response(m, k, alpha, beta, nfreq, fmin, fmax, fcn)
end program
```

The computed frequency response functions.

![](images/frf_proportional_example_1.png?raw=true)

## Nonlinear FRF Example
Computing the frequency response function for a nonlinear system is not as straight-forward.  A technique for capturing nonlinear behaviors, such as jump phenomenon, is to sweep through frequency, in both an ascending and a descending manner.  This example illustrates such a frequency sweeping using the famous Duffing equation as the model.

```math
\ddot{x} + \delta \dot{x} + \alpha x + \beta x^3 = \gamma \sin \omega t
```

The following module contains the equation.
```fortran
module duffing_ode_container
    use iso_fortran_env
    use dynamics
    implicit none

    ! Duffing Model Parameters
    real(real64), parameter :: alpha = 1.0d0
    real(real64), parameter :: beta = 4.0d-2
    real(real64), parameter :: delta = 1.0d-1
    real(real64), parameter :: gamma = 1.0d0

contains
    pure subroutine duffing_ode(freq, x, y, dydx, args)
        real(real64), intent(in) :: freq
            ! The excitation frequency
        real(real64), intent(in) :: x
            ! The independent variable.
        real(real64), intent(in), dimension(:) :: y
            ! An array of the N dependent variables.
        real(real64), intent(out), dimension(:) :: dydx
            ! An output array of length N where the derivatives are written.
        class(*), intent(inout), optional :: args
            ! An optional object for input/output of additional information.

        ! Variables
        real(real64) :: f

        ! Compute the harmonic forcing function
        f = gamma * sin(freq * x)

        ! Compute the derivatives
        dydx(1) = y(2)
        dydx(2) = f - delta * y(2) - alpha * y(1) - beta * y(1)**3
    end subroutine
end module
```
The calling program is as follows (plotting code ommitted).
```fortran
program example
    use iso_fortran_env
    use dynamics
    use duffing_ode_container
    implicit none

    ! Parameters
    real(real64), parameter :: f1 = 0.5d0
    real(real64), parameter :: f2 = 2.0d0
    integer(int32), parameter :: nfreq = 100
    
    ! Local Variables
    procedure(harmonic_ode), pointer :: fcn
    type(frf) :: solup, soldown

    ! Point to the ODE routine
    fcn => duffing_ode

    ! Perform the ascending sweep
    solup = frequency_sweep(fcn, nfreq, f1, f2, [0.0d0, 0.0d0])

    ! Perform the descending sweep
    soldown = frequency_sweep(fcn, nfreq, f2, f1, [0.0d0, 0.0d0])
end program
```
The computed frequency response functions, both ascending and descending, as compared with the analytical approximation.

![](images/frf_sweep_example_1.png?raw=true)

## Parameter Discovery (System Identification):
The following example illustrates how to estimate parameters of an ODE given an observed output to a known input.  This example illustrates how to find $\omega_{n}$ and $\zeta$ in the model of a single degree of freedom system.
```math
\ddot{x} + 2 \zeta \omega_{n} \dot{x} + \omega_{n}^{2} x = f(t)
```

```fortran
module equation_container
    use iso_fortran_env
    use dynamics
    implicit none

contains
    subroutine eom(t, x, dxdt, args)
        real(real64), intent(in) :: t               ! the current time value
        real(real64), intent(in) :: x(:)            ! the current state vector
        real(real64), intent(out) :: dxdt(:)        ! the derivatives
        class(*), intent(inout), optional :: args   ! model information

        ! Local Variables
        real(real64) :: zeta, wn, F

        ! Extract the model information
        select type (args)
        class is (model_information)
            wn = args%model(1)
            zeta = args%model(2)
            F = args%excitation%interpolate_value(t)
        end select

        ! The ODE:
        ! x" + 2 zeta wn x' + wn**2 x = F(t)
        dxdt(1) = x(2)
        dxdt(2) = F - wn * (2.0d0 * zeta * x(2) + wn * x(1))
    end subroutine
end module
```
The calling program is as follows (plotting code ommitted).
```fortran
program example
    use iso_fortran_env
    use equation_container
    use dynamics
    use diffeq
    use fstats
    implicit none

    ! Parameters
    real(real64), parameter :: fs = 1.0d3
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: zeta = 5.0d-2
    real(real64), parameter :: wn = 3.0d2
    real(real64), parameter :: sigma_pct = 1.0d-1
    real(real64), parameter :: amplitude = 1.0d4

    ! Local Variables
    integer(int32) :: i
    real(real64) :: dt, tmax, p(2), ic(2)
    type(dynamic_system_measurement) :: measurements(1)
    procedure(ode), pointer :: fcn
    type(ode_container) :: mdl
    type(runge_kutta_45) :: integrator
    type(model_information) :: info
    type(linear_interpolator), target :: interp
    real(real64), allocatable, dimension(:,:) :: sol
    type(iteration_controls) :: controls
    type(regression_statistics) :: stats(2)
    
    ! Generate an initial guess
    p = [2.5d2, 1.0d-1]

    ! Allocate memory for the measurement data we're trying to fit
    allocate( &
        measurements(1)%t(npts), &
        measurements(1)%output(npts), &
        measurements(1)%input(npts) &
    )

    ! Generate a time vector at which to sample the system.
    dt = 1.0d0 / fs
    tmax = dt * (npts - 1.0d0)
    measurements(1)%t = (/ (dt * i, i = 0, npts - 1) /)

    ! Define the forcing function at each time point
    measurements(1)%input = amplitude

    ! Generate the solution for the system
    ic = 0.0d0  ! zero-valued initial conditions
    call interp%initialize(measurements(1)%t, measurements(1)%input)
    info%model = [wn, zeta]
    info%excitation => interp
    mdl%fcn => eom
    call integrator%solve(mdl, measurements(1)%t, ic, args = info)
    sol = integrator%get_solution()
    measurements(1)%output = sol(:,2) + &
        box_muller_sample(0.0d0, 5.0d-3, npts) ! additional noise

    ! This is optional, but is illustrated here to show how to adjust solver
    ! tolerances
    call controls%set_to_default()
    controls%change_in_solution_tolerance = 1.0d-12
    controls%residual_tolerance = 1.0d-8

    ! Set up the problem and solve
    fcn => eom
    call siso_model_fit_least_squares(fcn, measurements, ic, p, &
        controls = controls, stats = stats)

    ! Compare the solution and the actual values
    print "(A)", "NATURAL FREQUENCY TERM:"
    print "(A,A,F8.3,A)", achar(9), "Actual: ", wn, " rad/s"
    print "(A,A,F8.3,A)", achar(9), "Computed: ", p(1), " rad/s"
    print "(A,A,F8.3,A)", achar(9), "Difference: ", p(1) - wn, " rad/s"
    print "(A,A,F8.3,A)", achar(9), "Std. Error: ", stats(1)%standard_error, " rad/s"
    print "(A,A,F8.3,A)", achar(9), "Conf. Int.: +/-", stats(1)%confidence_interval, " rad/s"
    print "(A,A,EN10.3)", achar(9), "P-Value: ", stats(1)%probability
    print "(A,A,EN12.3)", achar(9), "T-Statistic: ", stats(1)%t_statistic
    
    print "(A)", "DAMPING TERM:"
    print "(A,A,F6.3)", achar(9), "Actual: ", zeta
    print "(A,A,F6.3)", achar(9), "Computed: ", p(2)
    print "(A,A,F6.3)", achar(9), "Difference: ", p(2) - zeta
    print "(A,A,F6.3)", achar(9), "Std. Error: ", stats(2)%standard_error
    print "(A,A,F6.3)", achar(9), "Conf. Int.: +/-", stats(2)%confidence_interval
    print "(A,A,EN10.3)", achar(9), "P-Value: ", stats(2)%probability
    print "(A,A,EN12.3)", achar(9), "T-Statistic: ", stats(2)%t_statistic
end program
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

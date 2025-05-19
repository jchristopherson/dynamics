# dynamics
A library of routines used for the analysis of dynamic systems.

## Status
[![CMake](https://github.com/jchristopherson/dynamics/actions/workflows/cmake.yml/badge.svg)](https://github.com/jchristopherson/dynamics/actions/workflows/cmake.yml)
[![Actions Status](https://github.com/jchristopherson/dynamics/workflows/fpm/badge.svg)](https://github.com/jchristopherson/dynamics/actions)

## Documentation
The documentation can be found [here](https://jchristopherson.github.io/dynamics/).

## Capabilities
Here is a high-level list of the capabilities of this library.
- Compute linear frequency response functions for LTI systems.
- Perform modal analysis of an LTI system.
- Compute the frequency response of nonlinear systems in such a manner as to expose nonlinear behaviors such as jump phenomenon.
- Fit transfer functions to experimental data.
- Describe rigid body rotation and translation.
- Perform forward and inverse kinematic analysis for linkages.
- Analyze structural vibrations problems via linear 2D and 3D beam FEM.
- Determine properties of vibrating systems from experimental data such as resonant frequency, damping ratio, Q-factor, rise time, settling amplitudes, etc.
- Evaluate the step response behavior of SDOF systems.

## Kinematics Example
The following example illustrates the forward and inverse kinematic models of the illustrated 3R mechanism.  This example is Example 127 from Jazar's text "Theory of Applied Robotics, Kinematics, Dynamics, & Control."

![](images/3R%20Manipulator.PNG?raw=true)

The following module describes the geometry of the linkage.
```fortran
module linkage
    use iso_fortran_env
    use dynamics
    
    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Model Properties
    real(real64), parameter :: L1 = 1.5d0
    real(real64), parameter :: L2 = 2.0d1
    real(real64), parameter :: L3 = 1.0d1

    ! Denavit-Hartenberg Parameters
    real(real64) :: a(3) = [0.0d0, L2, L3]
    real(real64) :: alpha(3) = [0.5d0 * pi, 0.0d0, 0.0d0]
    real(real64) :: d(3) = [0.0d0, L1, 0.0d0]

contains
    subroutine kinematics_equations(jointvars, equations)
        ! The kinematics equations.
        real(real64), intent(in), dimension(:) :: jointvars
            ! The joint variables.
        real(real64), intent(out), dimension(:) :: equations
            ! The resulting kinematic equations.

        ! Local Variables
        real(real64) :: T(4, 4)
        
        ! Compute the forward kinematics problem
        T = dh_forward_kinematics(alpha, a, jointvars, d)

        ! Define the equations.
        ! 1. X position of the end effector
        ! 2. Y position of the end effector
        ! 3. Z position of the end effector
        ! 4. Orientation component of the end-effector
        ! 5. Orientation component of the end-effector
        ! 6. Orientation component of the end-effector
        equations(1:3) = T(1:3,4)
        equations(4) = T(1,1)
        equations(5) = T(2,2)
        equations(6) = T(3,1)
    end subroutine
end module
```
The kinematics code is as follows.
```fortran
program example
    use iso_fortran_env
    use dynamics
    use linkage
    implicit none

    ! Local Variables
    real(real64) :: theta(3), T(4, 4), qo(3), q(3), constraints(6)
    procedure(vecfcn), pointer :: mdl

    ! Define the Denavit-Hartenberg (DH) parameters
    call random_number(theta)   ! Randomly assign theta.  This is the joint variable

    ! Compute the forward kinematics problem
    T = dh_forward_kinematics(alpha, a, theta, d)

    ! -------------------------
    ! Solve the inverse problem.  Use the end-effector position and orientation
    ! computed by the forward kinematics process as a target for the inverse
    ! calculations.

    ! First define an initial guess
    qo = [0.0d0, 0.0d0, 0.0d0]

    ! Define the constraints for each kinematic equation
    constraints(1:3) = T(1:3,4)
    constraints(4) = T(1,1)
    constraints(5) = T(2,2)
    constraints(6) = T(3,1)

    ! Solve the model
    mdl => kinematics_equations
    q = solve_inverse_kinematics(mdl, qo, constraints)
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
The calling program is as follows.
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
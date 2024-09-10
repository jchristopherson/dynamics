# dynamics
A library of routines used for the analysis of dynamic systems.

## Status
[![CMake](https://github.com/jchristopherson/dynamics/actions/workflows/cmake.yml/badge.svg)](https://github.com/jchristopherson/dynamics/actions/workflows/cmake.yml)
[![Actions Status](https://github.com/jchristopherson/dynamics/workflows/fpm/badge.svg)](https://github.com/jchristopherson/dynamics/actions)

## Documentation
The documentation can be found [here](https://jchristopherson.github.io/dynamics/).

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
    pure subroutine duffing_ode(freq, x, y, dydx)
        real(real64), intent(in) :: freq
            ! The excitation frequency
        real(real64), intent(in) :: x
            ! The independent variable.
        real(real64), intent(in), dimension(:) :: y
            ! An array of the N dependent variables.
        real(real64), intent(out), dimension(:) :: dydx
            ! An output array of length N where the derivatives are written.

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

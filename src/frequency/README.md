## Frequency-Analysis Modules

* **`dynamics_modal_analysis`** solves the generalized eigenproblem for dense or CSR mass and stiffness matrices, returns modal frequencies and mode shapes, normalizes mode shapes, and computes proportional modal damping factors.
* **`dynamics_frequency_response`** computes single- and multiple-input frequency response functions, supports modal excitation, provides receptance/accelerance models and dynamic stiffness, and fits FRF data with nonlinear least squares.
* **`dynamics_frequency_sweep`** integrates ODE models under harmonic excitation across a frequency range, including ascending or descending sweeps for nonlinear response, and provides a linear chirp excitation function.

## Frequency Response Example
Consider the following 3 DOF system. The [`frf_example_1`](../../examples/frf_example_1.f90) example illustrates how to use this library to compute the frequency response functions for this system.

![](../../images/3%20DOF%20Schematic.PNG?raw=true)

The equations describing this system are as follows.

```math
\begin{bmatrix} m_1 & 0 & 0 \\ 0 & m_2 & 0 \\ 0 & 0 & m_3 \end{bmatrix} \begin{Bmatrix} \ddot{x}_1 \\ \ddot{x}_2 \\ \ddot{x}_3 \end{Bmatrix} + \begin{bmatrix} b_1 + b_2 & -b_2 & 0 \\ -b_2 & b_2 + b_3 & -b_3 \\ 0 & -b_3 & b_3 + b_4 \end{bmatrix} \begin{Bmatrix} \dot{x}_1 \\ \dot{x}_2 \\ \dot{x}_3 \end{Bmatrix} + \begin{bmatrix} k_1 + k_2 & -k_2 & 0 \\ -k_2 & k_2 + k_3 & -k_3 \\ 0 & -k_3 & k_3 + k_4 \end{bmatrix} \begin{Bmatrix} x_{1} \\ x_{2} \\ x_{3} \end{Bmatrix} = \begin{Bmatrix} F(t) \\ 0 \\ 0 \end{Bmatrix}
```

The essential excitation and solution setup is:

```fortran
real(real64), dimension(3,3) :: mass, stiffness, damp
type(frf) :: response
procedure(modal_excite), pointer :: excitation

mass = reshape([0.5d0, 0.0d0, 0.0d0, &
    0.0d0, 2.5d0, 0.0d0, &
    0.0d0, 0.0d0, 0.75d0], [3,3])
stiffness = reshape([15.0d6, -10.0d6, 0.0d0, &
    -10.0d6, 20.0d6, -10.0d6, &
    0.0d0, -10.0d6, 15.0d6], [3,3])
damp = reshape([150.0d0, -25.0d0, 0.0d0, &
    -25.0d0, 50.0d0, -25.0d0, &
    0.0d0, -25.0d0, 40.0d0], [3,3])

excitation => modal_frf_forcing_term
response = frequency_response(mass, damp, stiffness, &
    1000, 2.0d0*pi*10.0d0, 2.0d0*pi*1.0d3, excitation)

contains
subroutine modal_frf_forcing_term(freq, force, args)
    real(real64), intent(in) :: freq
    complex(real64), intent(out), dimension(:) :: force
    class(*), intent(in), optional :: args
    force = [(1.0d3, 0.0d0), (0.0d0, 0.0d0), (0.0d0, 0.0d0)]
end subroutine
```

The computed frequency response functions. The frequency values passed to `frequency_response` use the same units as the excitation callback; here the sweep bounds are angular frequencies in rad/s, corresponding to 10 Hz through 1 kHz.

![](../../images/frf_example_1.png?raw=true)

## Nonlinear FRF Example
Computing the frequency response function for a nonlinear system is not as straight-forward. A technique for capturing nonlinear behaviors, such as jump phenomenon, is to sweep through frequency, in both an ascending and a descending manner. The [`frf_sweep_example_1`](../../examples/frf_sweep_example_1.f90) example illustrates such a frequency sweep using the Duffing equation as the model.

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
    class(*), intent(in), optional :: args

    derivative(1) = state(2)
    derivative(2) = sin(frequency*t) - 0.1d0*state(2) - &
        state(1) - 0.04d0*state(1)**3
end subroutine
```
The computed frequency response functions, both ascending and descending, as compared with the analytical approximation. Different sweep directions can follow different steady-state response branches; this is the hysteresis/jump behavior the example is intended to show.

![](../../images/frf_sweep_example_1.png?raw=true)
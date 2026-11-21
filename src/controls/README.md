## Controls Module

The `dynamics_controls` module provides continuous-time linear control-system tools:

* `state_space` represents a system as $\dot{x}=Ax+Bu,\; y=Cx+Du$. It has constructors for general matrix systems, second-order mass-damping-stiffness systems, and PID-controlled plants. Its methods evaluate derivatives and outputs and compute poles, zeros, and transfer functions.
* `transfer_function` and `polynomial` represent and combine scalar transfer functions. Transfer functions can be evaluated in the Laplace domain or at a supplied angular frequency, and converted to controllable or observable canonical state-space form.
* `lti_solve` time-integrates a state-space system for a user-supplied excitation callback (`ss_excitation`). An optional ODE integrator and user data can also be supplied.

## PID Control Example

The [`pid_control_example_1`](../../examples/pid_control_example_1.f90) example closes a PID loop around a mass-spring-damper plant. It uses a square-wave reference, compares the closed-loop output to the command over 0.5 seconds, and evaluates the closed-loop transfer function over 1--1000 Hz.

The example's `model` module defines the plant parameters and the square-wave reference. The mass, damping, and stiffness values are $m=6.0\times10^{-4}$, $b=2\zeta\sqrt{km}$, and $k=2.0\times10^3$; $\zeta=0.3$. The command has amplitude $X=0.5$ and frequency $f=5$ Hz:

```fortran
subroutine command_signal(t, r, args)
    real(real64), intent(in) :: t
    real(real64), intent(out) :: r(:)
    class(*), intent(in), optional :: args

    r = sign(X, sin(2.0d0 * pi * f * t))
end subroutine
```

The core setup creates the second-order plant, wraps it in a closed-loop PID model, and integrates it from rest:

```fortran
use iso_fortran_env, only : real64
use dynamics, only : state_space, ss_excitation, lti_solve
use model, only : m, b, k, X, command_signal

real(real64), parameter :: Kp = 0.551048902d0
real(real64), parameter :: Ki = 1295.323281d0
real(real64), parameter :: Kd = 1.086432759d-3
real(real64), parameter :: tau = 1.369008811d-4
real(real64), parameter :: tmax = 0.5d0
type(state_space) :: plant, mdl
procedure(ss_excitation), pointer :: fcn
real(real64), dimension(4) :: ic
real(real64), allocatable, dimension(:,:) :: sol

plant = state_space(m, b, k)
plant%C = reshape([k * X, 0.0d0], [1, 2])
mdl = state_space(Kp, Ki, Kd, tau, plant)

fcn => command_signal
ic = [0.0d0, 0.0d0, 0.0d0, 0.0d0]
sol = lti_solve(mdl, fcn, [0.0d0, tmax], ic)
```

Here `Kp`, `Ki`, and `Kd` are the proportional, integral, and derivative gains. `tau` sets the derivative-filter time constant in the PID state-space model. The closed-loop model has four states: the plant position and velocity, the integral-error state, and the filtered derivative state. The example plots the output alongside the command; it also computes the model's poles and zeros.

The gains and filter time constant are chosen by matching the closed-loop characteristic polynomial to four real poles at $-1400$, $-1600$, $-2200$, and $-3200$ per second, using the example's scaled displacement output $y=kXq$. This improves nominal square-wave tracking, at the cost of a larger force spike at each command edge. The example does not model actuator limits or sensor noise.

![Closed-loop output and square-wave command](../../images/pid_control_example_1.png?raw=true)

The Bode plot is evaluated at angular frequencies $\omega=2\pi f$, where $f$ spans 1--1000 Hz. The code converts the transfer-function magnitude to dB with $20\log_{10}\lvert H\rvert$ and the phase to degrees.

![Closed-loop PID Bode magnitude and phase](../../images/pid_control_example_1_bode.png?raw=true)

Build with `BUILD_DYNAMICS_EXAMPLES=ON`, then build and run the `pid_control_example_1` target. The Bode image in the repository is named `pid_control_example_1_bode.png`.
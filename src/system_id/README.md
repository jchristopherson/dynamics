## Parameter Discovery (System Identification):
The [`siso_lsq_fit_example`](examples/siso_lsq_fit_example.f90) example illustrates how to estimate parameters of an ODE given an observed output to a known input. It finds $\omega_{n}$ and $\zeta$ in the model of a single degree of freedom system.
```math
\ddot{x} + 2 \zeta \omega_{n} \dot{x} + \omega_{n}^{2} x = f(t)
```

The `dynamics_system_id` module estimates the parameters of an ODE model from measured data. For each candidate parameter set, the fitting routine integrates the model, compares the simulated output with the measured output, and uses a Levenberg-Marquardt least-squares solver to update the parameters until the residuals are minimized.

### Functionality in `dynamics_system_id`

* `dynamic_system_measurement`: holds one data set through the allocatable arrays `t` (sample times), `input` (excitation), and `output` (measured response).
* `siso_model_fit_least_squares`: generic interface to the fitting routine. It has two specific versions:
    * Version 1 takes a single initial-condition vector that is shared by all data sets.
    * Version 2 takes an M-by-NEQN initial-condition matrix, with one row per data set. Each record is integrated from its own initial condition, and all residuals are combined into one objective.
* `model_information`: the container the fitting routine passes to the ODE routine through its `args` argument. It holds the current parameter estimate (`model`), an interpolator for the excitation (`excitation`), and any user data (`user_info`). Use `select type` to unpack it, as shown in the `eom` stub below.
* `regression_statistics`: per-parameter results that describe the quality of the fit, including the standard error, confidence interval, p-value, and t-statistic.
* `iteration_controls`, `lm_solver_options`: optional settings for convergence tolerances and solver behavior. Call `set_to_default()` on `iteration_controls` to start from defaults.
* `convergence_info`: optional output describing how the iteration converged.
* `iteration_update`: optional callback interface for monitoring each solver iteration.
* `constraint_equations`: interface for adding constraints beyond the differential equations, such as quasi-static measurements that help identify a stiffness term.
* `ode_integrator`: base type for the integrator used to solve the model. If none is supplied, the default is the Dormand-Prince Runge-Kutta method.

The remaining optional arguments of the fitting routine are summarized below.

| Argument | Purpose |
| --- | --- |
| `integrator` | Integrator used to solve the model equations. |
| `ind` | Index of the ODE solution component to fit (default `1`). |
| `minp`, `maxp` | Lower and upper bounds on the parameters. |
| `stats` | Receives one `regression_statistics` entry per parameter. |
| `alpha` | Significance level for confidence intervals (default `0.05`, a 95% interval). |
| `cov` | Receives the parameter covariance matrix. |
| `weights` | Weighting factors for each equation. |
| `xc`, `yc`, `constraints` | Additional constraint data and equations; must be used together. |
| `args` | User data made available to the model through `user_info`. |

### Example

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
    class(*), intent(in), optional :: args
    ! Extract wn, zeta, and forcing from args; then evaluate the SDOF model.
end subroutine
```
The snippet shows the key steps of the full [example](examples/siso_lsq_fit_example.f90):

* **Measurements:** The example generates its own data by integrating the model with a 1 kHz sample rate, 1000 points (a 1 second record), and a constant forcing of `1.0d4`. It then adds Gaussian noise with a standard deviation of `5.0d-3` to the output. In practice, `t`, `input`, and `output` would be loaded from recorded data.
* **Initial guess:** `parameters` holds the starting estimate of `[wn, zeta]`. The solver updates this array in place, so on return it contains the fitted values. The estimate does not need to be exact, but it should be in a reasonable range.
* **Initial conditions:** `initial_state` is the state vector at the first time point. Because a single vector is passed, this is the Version 1 form of the fitting routine.
* **Model routine:** `model` is a procedure pointer to the ODE routine. On each iteration, the fitting routine fills `model_information` with the current parameters and an interpolator of the measured input, so `eom` can read `wn`, `zeta`, and the forcing value at time `t`. The `args` argument of the fitting routine is not needed for this; it only supplies `user_info`.
* **Statistics and controls:** `stats` receives one entry per parameter. The `controls` object is optional. The example tightens `change_in_solution_tolerance` and `residual_tolerance` to show how to adjust convergence. The defaults are used otherwise.
* **Verification:** After the fit, the example clears the integrator buffer and re-solves the model with the fitted parameters. It then plots the measured data against the fitted response.

The results are as follows. In the output, "Difference" is the computed value minus the actual value. The confidence interval is the 95% interval for each estimate. A small standard error and a large t-statistic indicate that the data constrain the parameter well, while a small p-value indicates the estimate is statistically distinguishable from zero. The differences shown here are consistent with the added noise.
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
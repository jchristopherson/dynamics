#include <stdio.h>
#include <math.h>
#include "dynamics.h"

void excitation(int n, double t, double *u);

int main()
{
    // Local Variables
    const double m = 1.0;
    const double b = 5.0;
    const double k = 2.0e3;
    const double tmax = 1.0;
    const int n = 100;
    int i;
    double ic[2], t[n], y[n], dt;
    c_state_space_model mdl;

    // Create a new model
    c_create_state_space_model(1, 1, &m, 1, &b, 1, &k, 1, &mdl);

    // Define the solution points
    dt = tmax / (n - 1.0);
    for (i = 0; i < n; ++i) t[i] = i * dt;

    // Define the initial conditions
    ic[0] = 0.0;
    ic[1] = 0.0;
    
    // Compute the solution
    c_lti_solve(&mdl, excitation, n, t, 2, ic, DYN_RUNGE_KUTTA_45, 1, y, n);

    // Print out the solution
    for (i = 0; i < n; ++i)
    {
        printf("%0.3f\t%0.3f\n", t[i], y[i]);
    }

    // End
    c_free_state_space_model(&mdl);
    return 0;
}

void excitation(int n, double t, double *u)
{
    u[0] = 1.0e3 * sin(10.0 * t);
}

// 3. 2-point boundary value problem by shooting.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void derivs (double x, double y[], double dydx[]);
void rk4 (double x, double y[], double dydx[], double h, double y_out[], void derivs (double x, double y[], double dydx[]));


int main() {
    FILE *filename;
    char name_file[150];
    double t = 0, h = 0, t_min = 0, t_max = 0, delta = 0, delta_min = 0, delta_max = 0, err = 10;
    double y[2], dydx[2], y_out[2];
    int i = 0, t_step = 0;
    
    // Number of steps.
    t_min = 0.;
    t_max = 1.;
    t_step = 1000;
    h = (t_max - t_min) / t_step;
    
    // Initial conditions.
    y[0] = 0;
    delta_min = 0.;
    delta_max = 5.;
    t = t_min;
    
    // Iterate until error smaller than some number.
    while (err > 1e-8) {
        // Use bisection to iterate the root of delta.
        delta = (delta_max + delta_min) / 2.;
        y[1] = delta;
        
        t = t_min;
        // Shooting.
        while (t <= t_max) {
            derivs (t, y, dydx);
            rk4 (t, y, dydx, h, y_out, derivs);
            t += h;
            y[0] = y_out[0];
            y[1] = y_out[1];
        }
        // Check the relative error of the trial at t = 1.
        err = fabs((y[0] - 1.) / 1.);
        
        // For next iteration.
        if (y[0] < 1)
            delta_min = delta;
        else
            delta_max = delta;
    }
    
    // After the iteration, the lower and upper limits of delta are fixed.
    t = t_min;
    y[0] = 0;
    y[1] = (delta_max + delta_min) / 2.;
    
    strcpy (name_file, "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps6/data/q3/iterated.dat");
    
    filename = fopen (name_file, "w");
    
    fprintf(filename, "%f\t%f\t%f\n", t, y[0], y[1]);
    
    // Print out iterated solution.
    while (t < t_max) {
        derivs (t, y, dydx);
        rk4 (t, y, dydx, h, y_out, derivs);
        // Going to the next step.
        t += h;
        y[0] = y_out[0];
        y[1] = y_out[1];
        fprintf(filename, "%f\t%f\t%f\n", t, y_out[0], y_out[1]);
    }
    
    fclose (filename);
    
    return 0;
}


void derivs (double x, double y[], double dydx[]) {
    // A function to store up the drivatives at each step.
    dydx[0] = y[1];
    dydx[1] = -(y[0] + 1) * M_PI * M_PI / 4.;
}

void rk4 (double x, double y[], double dydx[], double h, double y_out[], void derivs (double x, double y[], double dydx[])) {
    // 4th order Runge-Kutta.
    int i = 0, n = 2;
    double xhh = 0, xh = 0, hh = 0;
    
    // Half step.
    hh = h * 0.5;
    xh = x + h;
    xhh = x + hh;
    
    double *dym = malloc (n * sizeof (double));
    double *dyt = malloc (n * sizeof (double));
    double *dyf = malloc (n * sizeof (double));
    double *yt = malloc (n * sizeof (double));
    
    for (i = 0; i < n; i++) {
        yt[i] = y[i] + hh * dydx[i];    // First step, dydx[i] is k1.
    }
    derivs (xhh, yt, dyt);    // k2
    
    for (i = 0; i < n; i++) {
        yt[i] = y[i] + hh * dyt[i];    // Second step.
    }
    derivs (xhh, yt, dym);    // k3
    
    for (i = 0; i < n; i++) {
        yt[i] = y[i] + h * dym[i];
    }
    derivs (xh, yt, dyf);    // k4

    for (i = 0; i < n; i++) {
        y_out[i] = y[i] + (h / 6.) * (dydx[i] + 2 * (dyt[i] + dym[i]) + dyf[i]);
    }
}
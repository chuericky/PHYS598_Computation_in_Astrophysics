// 1. Nonlinear dynamics and chaos.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void derivs (double x, double y[], double dydx[], double q, double b, double w0);
void rk4 (double x, double y[], double dydx[], double h, double y_out[], void derivs (double x, double y[], double dydx[], double q, double b, double w0), double q, double b, double w0);

int main() {
    FILE *filename;
    char name_file[150];
    double t = 0, h = 0, t_min = 0, t_max = 0, q = 0, b = 0, w0 = 0;
    double y[2], dydx[2], y_out[2];
    int i = 0, t_step = 0;
    
    
    // Initial conditions.
    y[0] = 0;
    y[1] = 2;
    t = t_min;
    
    // Parameters.
    q = 0.5;
    b = 1.15;    //1.15
    w0 = 2./3;
    derivs (t, y, dydx, q, b, w0);
    
    // Number of steps.
    t_min = 0.;
    t_max = 5000.;
    t_step = 10000;
    h = (t_max - t_min) / t_step;
    
    strcpy (name_file, "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps6/data/q1/4.dat");
    
    filename = fopen (name_file, "w");
    
    while (t < t_max) {
        rk4 (t, y, dydx, h, y_out, derivs, q, b, w0);
        fprintf(filename, "%f\t%f\t%f\n", t, y_out[0], y_out[1]);
        // Going to the next step.
        t += h;
        y[0] = y_out[0];
        y[1] = y_out[1];
        derivs (t, y, dydx, q, b, w0);
    }
    
    fclose (filename);
    
    return 0;
}

void derivs (double x, double y[], double dydx[], double q, double b, double w0) {
    // A function to store up the drivatives at each step.
    dydx[0] = y[1];
    dydx[1] = -q * y[1] - sin(y[0]) + b * cos(w0 * x);
}

void rk4 (double x, double y[], double dydx[], double h, double y_out[], void derivs (double x, double y[], double dydx[], double q, double b, double w0), double q, double b, double w0) {
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
    derivs (xhh, yt, dyt, q, b, w0);    // k2
    
    for (i = 0; i < n; i++) {
        yt[i] = y[i] + hh * dyt[i];    // Second step.
    }
    derivs (xhh, yt, dym, q, b, w0);    // k3
    
    for (i = 0; i < n; i++) {
        yt[i] = y[i] + h * dym[i];
    }
    derivs (xh, yt, dyf, q, b, w0);    // k4

    for (i = 0; i < n; i++) {
        y_out[i] = y[i] + (h / 6.) * (dydx[i] + 2 * (dyt[i] + dym[i]) + dyf[i]);
    }
}
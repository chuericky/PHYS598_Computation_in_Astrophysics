// 2. Stiff system of equations, using RK scheme.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void derivs (double x, double y[], double dydx[]);
void rk4 (double x, double y[], double dydx[], double h, double y_out[], void derivs (double x, double y[], double dydx[]));

int main() {
    FILE *filename;
    char name_file[150];
    double t = 0, h = 0, t_min = 0, t_max = 0;
    double y[3], dydx[3], y_out[3];
    int i = 0, t_step = 0;
    
    // Initial conditions.
    y[0] = 1;
    y[1] = 1;
    y[2] = 0;
    t = t_min;
    derivs (t, y, dydx);
    
    // Number of steps.
    t_min = 0.;
    t_max = 100.;
    t_step = 1000;
    h = (t_max - t_min) / t_step;
    
    strcpy (name_file, "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps6/ps6_tex/q2/rk4.dat");
    
    filename = fopen (name_file, "w");
    
    while (t < t_max) {
        rk4 (t, y, dydx, h, y_out, derivs);
        fprintf(filename, "%f\t%f\t%f\t%f\n", t, y_out[0], y_out[1], y_out[2]);
        // Going to the next step.
        t += h;
        y[0] = y_out[0];
        y[1] = y_out[1];
        y[2] = y_out[2];
        derivs (t, y, dydx);
    }
    
    fclose (filename);
    
    return 0;
}


void derivs (double x, double y[], double dydx[]) {
    // A function to store up the drivatives at each step.
    dydx[0] = -0.013 * y[0] - 1000 * y[0] * y[2];
    dydx[1] = -2500 * y[1] * y[2];
    dydx[2] = -0.013 * y[0] - 1000 * y[0] * y[2] - 2500 * y[1] * y[2];
}

void rk4 (double x, double y[], double dydx[], double h, double y_out[], void derivs (double x, double y[], double dydx[])) {
    // 4th order Runge-Kutta.
    int i = 0, n = 3;
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
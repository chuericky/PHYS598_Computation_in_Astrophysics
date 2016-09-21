// 5. 1-D Shrodinger Equation

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


void derivs (double x, double y[], double dydx[], double E);
void rk4 (double x, double y[], double dydx[], double h, double y_out[], void derivs (double x, double y[], double dydx[], double E), double E);
double trapezoidal (int N, double a, double b, double u[]);


int main() {
    FILE *filename;
    char name_file[150];
    double t = 0, h = 0, t_min = 0, t_max = 0, t_match = 0, E = 0, E_min = 0, E_max = 0, V = 0, left_t = 0, right_t = 0, err = 10, accept = 1e-5, norm = 0;
    double y[2], dydx[2], y_out[2];
    int i = 0, k = 0, t_step = 0;
    
    // Initial conditions.
    t = t_min;
    E_min = 2;    // Strict lower bound of eigen-energy.
    E_max = 4;
    t_step = 1000;
    E = E_min;
    
    double **y_left = malloc (2 * sizeof (double *));
    double **y_right = malloc (2 * sizeof (double *));
    double **y_all = malloc (2 * sizeof (double *));
    double **y_all_sq = malloc (2 * sizeof (double *));
    
    for (i = 0; i < 2; i++) {
        y_left[i] = malloc((t_step + 1) * sizeof(double));
        y_right[i] = malloc((t_step + 1) * sizeof(double));
        y_all[i] = malloc((2 * t_step + 1) * sizeof(double));
        y_all_sq[i] = malloc((2 * t_step + 1) * sizeof(double));
        
        for (k = 0; k < t_step + 1; k++) {
            y_left[i][k] = 0;
            y_right[i][k] = 0;
        }
        for (k = 0; k < 2 * t_step + 1; k++) {
            y_all[i][k] = 0;
            y_all_sq[i][k] = 0;
        }
    }
    
    strcpy (name_file, "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps6/data/q5/E3.dat");
    
    filename = fopen (name_file, "w");
    
    // Iterate until error goes below the tolerance.
    while (err > accept) {
        
        // Integration from the left to the right.
        // Initial conditions.
        y_left[0][0] = 0;
        y_left[1][0] = 1;    // Without loss of generality, take u_left'(0) = 1.
        y[0] = y_left[0][0];
        y[1] = y_left[1][0];
    
        // At the end-points of integration.
        t_match = acosh(sqrt(6 / (3 - E)));
        t_min = t_match - fabs(t_match) * (t_step / 100.0);
        t_max = t_match;
        h = (t_max - t_min) / t_step;
    
        t = t_min;
    
        // Do the integration.
        for (i = 1; i <= t_step; i++) {
            derivs (t, y, dydx, E);
            rk4 (t, y, dydx, h, y_out, derivs, E);
            t = t_min + i * h;
            // Store up the values.
            y_left[0][i] = y_out[0];
            y_left[1][i] = y_out[1];
            y[0] = y_out[0];
            y[1] = y_out[1];
        }
        left_t = y_out[1] / y_out[0];
    
    
        // Integration from the right to the left.
        for (k = 0; k < 2; k++) {
            // Initial conditions.
            y[0] = y_right[0][t_step] = 0;
        
            if (k == 0) {
                y[1] = y_right[1][t_step] = 1;
            }
            else {
                y[1] = y_right[1][t_step] = -1;
            }
        
            // At the end-points of integration.
            t_match = acosh(sqrt(6 / (3 - E)));
            t_min = t_match;
            t_max = t_match + fabs(t_match) * (t_step / 100.0);
            h = -(t_max - t_min) / t_step;
        
            t = t_max;
            
            // Integrate.
            for (i = 1; i <= t_step; i++) {
                derivs (t, y, dydx, E);
                rk4 (t, y, dydx, h, y_out, derivs, E);
                t = t_max + i * h;
                // Store up the values.
                y_right[0][t_step - i] = y_out[0];
                y_right[1][t_step - i] = y_out[1];
                y[0] = y_out[0];
                y[1] = y_out[1];
            }
            right_t = y_out[1] / y_out[0];
            
            err = fabs(left_t - right_t);
            
            if (err < accept) {
                break;
            }
        }
    
        
        if (err < accept) {
            for (i = 0; i <= 2 * t_step; i++) {
                if (i < t_step) {
                    // Store up the arrays on the left.
                    y_all[0][i] = y_left[0][i];
                    y_all[1][i] = y_left[1][i];
                }
                else {
                    y_all[0][i] = y_right[0][i - t_step] * y_left[0][t_step] / y_right[0][0];
                    y_all[1][i] = y_right[1][i - t_step] * y_left[0][t_step] / y_right[0][0];
                }
            }
            
            t_min = t_match - fabs(t_match) * (t_step / 100.0);
            h = (t_max - t_min) / (2 * t_step);
            
            for (i = 0; i <= 2 * t_step; i++) {
                y_all_sq[0][i] = pow(y_all[0][i], 2);
                y_all_sq[1][i] = pow(y_all[1][i], 2);
            }
            
            // Normalization factor.
            norm = pow((trapezoidal (2 * t_step, t_min, t_min + 2 * t_step * h, y_all_sq[0])), 0.5);
            
            for (i = 0; i <= 2 * t_step; i++) {
                t = t_min + i * h;
                V = 6 * (0.5 - 1 / pow(cosh(t), 2));
                fprintf(filename, "%f\t%f\t%f\t%f\t%f\n", t, y_all[0][i] / norm, y_all[1][i] / norm, E, V);
            }
        }
    
        else {
            E += (E_max - E_min) / 100.0;
            if (E > E_max) {
                printf("E > E_max!\n");
                exit(1);
            }
        }
    }
    
    free (y_left);
    free (y_right);
    free (y_all);
    fclose (filename);
    
    return 0;
}


void derivs (double x, double y[], double dydx[], double E) {
    // A function to store up the drivatives at each step.
    dydx[0] = y[1];
    dydx[1] = (12 * (0.5 - 1/pow(cosh(x), 2)) - 2 * E) * y[0];
}

void rk4 (double x, double y[], double dydx[], double h, double y_out[], void derivs (double x, double y[], double dydx[], double E), double E) {
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
    derivs (xhh, yt, dyt, E);    // k2
    
    for (i = 0; i < n; i++) {
        yt[i] = y[i] + hh * dyt[i];    // Second step.
    }
    derivs (xhh, yt, dym, E);    // k3
    
    for (i = 0; i < n; i++) {
        yt[i] = y[i] + h * dym[i];
    }
    derivs (xh, yt, dyf, E);    // k4

    for (i = 0; i < n; i++) {
        y_out[i] = y[i] + (h / 6.) * (dydx[i] + 2 * (dyt[i] + dym[i]) + dyf[i]);
    }
}

// Function for Trapezoidal's Rule.
double trapezoidal (int N, double a, double b, double u[]) {
    
    double h = 0, sum = 0;
    int i = 0;
    
    h = (b - a) / N;
    
    // Trapezoidal's Rule.
    if (N == 1) {
        sum += h * (u[0] + u[N]) / 2.;
    }
    
    // Extended Trapezoidal's Rule. (N > 1)
    else {
        sum += h * (u[0] + u[N]) / 2.;
        for (i = 0; i < N - 1; i++)
            sum += h * u[i + 1];
    }
    
    return sum;
}
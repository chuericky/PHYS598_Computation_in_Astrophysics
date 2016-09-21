// Q1.b.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "math.h"

int main() {
    int N_time = 1001;
    int N_space = 1001;
    int i = 0, j = 0;
    
    float x0 = -10.0, x1 = 10.0;
    float *a, *b, *c, *r, *u;
    
    double dx = (x1 - x0) / (N_space - 1.);     // Step in x.
    double dt = 0.025 * dx * dx;                  // Time step.
    double alpha = dt / (dx * dx);
    
    a = vector(0, N_space);
    b = vector(0, N_space);
    c = vector(0, N_space);
    r = vector(0, N_space);
    u = vector(0, N_space);
    
    // Initial condition.
    for (i = 0; i < N_space; i++) {
        u[i] = exp(-pow((x0 + i * dx), 2));
    }
        
    // The tridiagonal matrix.
    for (i = 0; i < N_space; i++) {
        if (i == 0 || i == N_space - 1) {
            a[i] = 0;
            b[i] = 0;
            c[i] = 0;
        }
        else if (i == 1) {
            a[i] = 0;
            b[i] = 1 + 2 * alpha;
            c[i] = -alpha;
        }
        else if (i == N_space - 2) {
            a[i] = -alpha;
            b[i] = 1 + 2 * alpha;
            c[i] = 0;
        }
        else {
            a[i] = -alpha;
            b[i] = 1 + 2 * alpha;
            c[i] = -alpha;
        }
    }
    
    for (j = 1; j < N_time; j++) {
        u[0] = 0;
        u[N_space - 1] = 0;
        
        for (i = 0; i < N_space; i++) {
            if (i == 0 || i == N_space - 1) r[i] = 0;
            else    r[i] = u[i];
        }
        
        tridag(a, b, c, r, u, N_space - 2);
        
        if (j == 900) {
            for (i = 0; i < N_space; i++) {
                printf("%d\t%d\t%f\t%f\n", j, i, x0 + i * dx, u[i]);
            }
        }
    }
    
    free_vector(a, 0, N_space);
    free_vector(b, 0, N_space);
    free_vector(c, 0, N_space);
    free_vector(r, 0, N_space);
    free_vector(u, 0, N_space);
    
    return 0;
}
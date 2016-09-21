// Q.1. Eigenvalues of a tridiagonal matrix.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "math.h"

int main() {
    int i = 0, j = 0;
    int N = 32;  // Size of the matrix
    double *d, *e, *f, **a;
    
    // The tridiagonal N * N matrix.
    double **C = malloc (N * sizeof (double *));
    
    for (i = 0; i < N; i++) {
        C[i] = malloc (N * sizeof (double));
    }
    
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            // The tridiagonal matrix.
            if (i == j)                 C[i][j] = 2;
            else if (fabs(i - j) == 1)  C[i][j] = -1;
            else                        C[i][j] = 0;
        }
    }
    
    d = dvector(1, N);
    e = dvector(1, N);
    f = dvector(1, N);
    a = dmatrix(1, N, 1, N);
    
    for (i = 1; i <= N; i++) {
        for (j = 1; j <= N; j++) {
            a[i][j] = C[i - 1][j - 1];
        }
    }
    
    // Transform the matrix with the orthogonal matrix Q.
    tred2(a, N, d, e);
    // QL algorithm.
    tqli(d, e, N, a);
    
    double *ana_ev = malloc (N * sizeof (double));      // Analytic eigenvalues.
    double *num_ev = malloc (N * sizeof (double));      // Numerical eigenvalues.
    
    for (i = 0; i < N; i++) {
        ana_ev[i] = 4 * pow(sin((i + 1) * M_PI / (2 * (N + 1))), 2);
        num_ev[i] = d[i + 1];
    }
    
    printf("Lambda_n\tNumerical\tAnalytic\tRelative Error\n");
    for (i = 0; i < N; i++) {
        printf("%d\t%.*lf\t%.*lf\t%lE\n", i + 1, 14, num_ev[i], 14, ana_ev[i], (num_ev[i] - ana_ev[i])/ana_ev[i]);
    }
    
    free_dmatrix (a, 1, N, 1, N);
    free_dvector (f, 1, N);
    free_dvector (d, 1, N);
    free_dvector (e, 1, N);
    free (C);
    free (ana_ev);
    free (num_ev);
    
    return 0;
}
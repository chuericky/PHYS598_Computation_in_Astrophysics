// Q.2. Eigenvalues of a tridiagonal matrix.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "math.h"

int main() {
    int i = 0, j = 0, k = 0, count = 0;
    int N = 32;  // Size of the matrix
    double *d, *e, *f, **a, normal = 0;
    
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
    
    // Eigenvectors
    for (i = 1; i <= N; i++) {
        for (j = 1; j <= N; j++) {
            f[j] = 0.;
            for (k = 1; k <= N; k++) {
                f[j] += (C[j - 1][k - 1] * a[k][i]);
            }
        }
        printf("\nNumerical Eigenvalue %d, omega^2 = %.*lf, i.e. omega = %.*lf\n\n", i, 14, d[i], 14, pow(d[i], 0.5));
        printf("%s\n", "Numerical Eigenvector");
        count = 0;
        for (j = 1; j <= N; j++) {
            printf("%.*lE\n", 14, a[j][i]);
            count++;
            /*if (count % 5 == 4) {
                count = 0;
                printf("\n");
            }*/
        }
        printf("\n");
    }
    
    // Normalization and the eigenvectors of the analytical solution.
    double *norm = malloc (N * sizeof (double));
    double **ana_vec = malloc (N * sizeof (double *));
    
    for (i = 0; i < N; i++) {
        norm[i] = 0;
        ana_vec[i] = malloc (N * sizeof (double));
        for (j = 0; j < N; j++) {
            ana_vec[i][j] = 0;
            norm[i] += pow(sin((j + 1) * (i + 1) * M_PI / (N + 1)), 2);
        }
        normal = -pow(1/norm[i], 0.5);
        for (j = 0; j < N; j++) {
            ana_vec[i][j] = normal * sin((j + 1) * (i + 1) * M_PI / (N + 1));
        }
        
        ana_ev[i] = 4 * pow(sin((i + 1) * M_PI / (2 * (N + 1))), 2);
    }
    
    for (i = 0; i < N; i++) {
        printf("\nAnalytical Eigenvalue %d, omega^2 = %.*lf, i.e. omega = %.*lf\n\n", i + 1, 14, ana_ev[i], 14, pow(ana_ev[i], 0.5));
        printf("%s\n", "Analytical Eigenvector");
        count = 0;
        for (j = 0; j < N; j++) {
            printf("%.*lE\n", 14, ana_vec[i][j]);
            count++;
            /*if (count % 5 == 4) {
                count = 0;
                printf("\n");
            }*/
        }
    }
    printf("\n");
    
    free_dmatrix (a, 1, N, 1, N);
    free_dvector (f, 1, N);
    free_dvector (d, 1, N);
    free_dvector (e, 1, N);
    free (C);
    free (ana_ev);
    free (norm);
    free (ana_vec);
    
    return 0;
}
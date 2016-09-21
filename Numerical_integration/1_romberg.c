// Q.1. Romberg's Method.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
    int i = 0, j = 0, N = 0, M = 0;
    double h = 0, sum = 0;
    
    do {
        // Input the number of steps. N >= 2, and N % 2 == 0.
        printf("Enter N then M (N >= M). R(N, M): ");
        scanf("%d %d", &N, &M);
    }while (N < M);
        
    M = M + 1;
    N = N + 1;
    
    h = (M_PI / 2.) - 0;
    
    // Dynamic array for Romberg's method.
    double **R = malloc(N * sizeof(double *));
    
    for (i = 0; i < N; i++) {
        R[i] = malloc(M * sizeof(double));
    }
    
    for (i = 0; i < N; i++) {
        for (j = 0; j < M; j++) {
            R[i][j] = 0;
        }
    }
    
    R[0][0] = h * (sin(M_PI/2.) + sin(0)) / 2.;
    printf(" R[0][0] = %.7f\n", R[0][0]);
    
    if (M == 1 & N == 1) {
        printf(" Iterated solution = %f.\n", R[0][0]);
    }
    
    else {
        // Romberg's method.
        for (i = 1; i < N; i++) {
            h *= 0.5;
            sum = 0;
        
            for (j = 1; j <= pow(2, (i - 1)); j++) {
                sum += sin((2 * j - 1) * h);
            }
            R[i][0] = 0.5 * R[i - 1][0] + sum * h;  // Trapezoidal rule with 2^N + 1 points.
            printf(" R[%d][0] = %.7f\n", i, R[i][0]);
            if (M != 1) {
                for (j = 1; j < M; j++) {
                    R[i][j] = R[i][j - 1] + (R[i][j - 1] - R[i - 1][j - 1]) / (pow(4, j) - 1);
                    printf(" R[%d][%d] = %.7f\n", i, j, R[i][j]);
                }
            }
        }
    }
    
    printf("\nIterated solution, R[%d][%d] = %.7lf.\n", N - 1, M - 1, R[N - 1][M - 1]);
    
    free(R);
    
    return 0;
}
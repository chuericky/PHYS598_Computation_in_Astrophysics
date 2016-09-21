// Q.1. Trapezoidal Rule.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
    int i = 0, N = 0, k = 0;
    double h = 0, sum = 0;
    
    // Input the number of steps. N >= 1.
    do {
        printf("Enter number of steps N, where N >= 1: ");
        scanf("%d", &N);
    }while (N < 1);
    
    h = (M_PI / 2.) / N;
    
    // Trapezoidal Rule.
    if (N == 1) {
        sum += h * (sin(0) / 2. + sin(M_PI/2.) / 2.);
    }
    // Extended Trapezoidal Rule. (N > 1)
    else {
        // Create an array
        double *sin_x = malloc((N - 1) * sizeof(double));
        for (i = 0; i < N - 1; i++) {
            sin_x[i] = sin(h * (i + 1));
        }
        
        sum += h * (sin(0) / 2. + sin(M_PI/2.) / 2.);
        for (i = 0; i < N - 1; i++) {
            sum += h * sin_x[i];
        }
        
        free(sin_x);
    }
    
    printf("N = %d\th = %lf\tI = %.7lf\n", N, h, sum);
    
    return 0;
}
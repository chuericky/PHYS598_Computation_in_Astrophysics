// Q.1. Simpson's Rule.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {
    int i = 0, N = 0;
    double h = 0, sum = 0;
    
    // Input the number of steps. N >= 2, and N % 2 == 0.
    do {
        printf("Enter number of steps N, where N >= 2, and N must be an even number: ");
        scanf("%d", &N);
    }while (N % 2 == 1);
    
    h = (M_PI / 2.) / N;
    
    // Simpson's Rule.
    if (N == 2) {
        sum += h * (sin(0) / 3. + 4 * sin(M_PI/4.) / 3. + sin(M_PI/2.) / 3.);
    }
    // Extended Simpson's Rule. (N > 2)
    else {
        // Create an array for storing sin(x).
        double *sin_x = malloc((N - 1) * sizeof(double));
        
        sum += h * (sin(0) / 3. + sin(M_PI/2.) / 3.);
        for (i = 0; i < N - 1; i++) {
            if (i % 2 == 0) {
                sin_x[i] = 4 * sin(h * (i + 1)) / 3.;
            }
            else {
                sin_x[i] = 2 * sin(h * (i + 1)) / 3.;
            }
        }
        
        for (i = 0; i < N - 1; i++) {
            sum += h * sin_x[i];
        }
        
        free(sin_x);
    }
    
    printf("N = %d\th = %lf\tI = %.7lf\n", N, h, sum);
    
    return 0;
}
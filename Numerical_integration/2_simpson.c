// Q.2. Simpson's Rule, without change of variable.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double f(double x);

int main() {
    int i = 0, N = 0;
    double h = 0, sum = 0;
    
    // Input the number of steps. N >= 2, and N % 2 == 0.
    do {
        printf("Enter number of steps N, where N >= 2, and N must be an even number: ");
        scanf("%d", &N);
    }while (N % 2 == 1);
    
    h = 5. / N;
    
    // Simpson's Rule.
    if (N == 2) {
        sum += h * (f(0) / 3. + 4 * f(2.5) / 3. + f(5) / 3.);
    }
    // Extended Simpson's Rule. (N > 2)
    else {
        // Create an array for storing sin(x).
        double *f_x = malloc((N - 1) * sizeof(double));
        
        sum += h * (f(0) / 3. + f(5) / 3.);
        for (i = 0; i < N - 1; i++) {
            if (i % 2 == 0) {
                f_x[i] = 4 * f(h * (i + 1)) / 3.;
            }
            else {
                f_x[i] = 2 * f(h * (i + 1)) / 3.;
            }
        }
        
        for (i = 0; i < N - 1; i++) {
            sum += h * f_x[i];
        }
        
        free(f_x);
    }
    
    printf("N = %d\th = %lf\tI = %.7lf\n", N, h, sum);
    
    return 0;
}


// The function 1/(|1.5 - x|)^0.5.
double f(double x) {
    return 1./sqrt(fabs(1.5 - x));
}
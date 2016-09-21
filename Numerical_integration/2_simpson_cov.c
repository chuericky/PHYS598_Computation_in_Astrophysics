// Q.2. Simpson's Rule, with change of variable.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int bool_dowhile(int N, int M);
double f();
double simpson(int N, double a, double b);

int main() {
    int i = 0, N = 0, M = 0;
    double h = 0, sum = 0;
    
    // There are 2 integrals. Input the number of steps for each integral. N, M >= 2, and N % 2 == 0 and M % 2 == 0.
    do {
        printf("Enter number of steps N and M, where N, M >= 2, and N and M must be even numbers: ");
        scanf("%d %d", &N, &M);
    }while (bool_dowhile(N, M) == 0);
    
    // The integration.
    sum = 2 * (simpson(N, 0, pow(1.5, 0.5)) + simpson(M, 0, pow(3.5, 0.5)));
    
    printf("N = %d\tM = %d\tI = %.7lf\n", N, M, sum);
    
    return 0;
}


// Boolean test for the do-while loop. This is to make sure N, M are multiples of 2.
int bool_dowhile(int N, int M) {
    if ((M % 2 == 0) && (N % 2 == 0) && (M != 0) && (N != 0))   return 1;
    else    return 0;
}

// Function in the integrand.
double f(double x) {
    return 1;
}

// Function for Simpson's Rule.
double simpson(int N, double a, double b) {
    
    double h = 0, sum = 0;
    int i = 0;
    
    h = (b - a) / N;
    
    // Simpson's Rule.
    if (N == 2) {
        sum += h * (f(a) / 3. + 4 * f(0.5 * (a + b)) / 3. + f(b) / 3.);
    }
    
    // Extended Simpson's Rule. (N > 2)
    else {
        // Create an array for storing sin(x).
        double *f_x = malloc((N - 1) * sizeof(double));
    
        sum += h * (f(a) / 3. + f(b) / 3.);
        for (i = 0; i < N - 1; i++) {
            if (i % 2 == 0) {
                f_x[i] = 4 * f(a + h * (i + 1)) / 3.;
            }
            else {
                f_x[i] = 2 * f(a + h * (i + 1)) / 3.;
            }
        }
    
        for (i = 0; i < N - 1; i++) {
            sum += h * f_x[i];
        }
        
        free(f_x);
    }
    
    return sum;
}
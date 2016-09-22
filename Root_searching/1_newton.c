// Q.1. Newton's method on finding the real roots of x^3 + 1.5x^2 - 5.75x + 4.37 = 0.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double f(double x);
double df(double x);

int main () {
    int precision = 0, i = 0;
    double x0 = 0, x1 = 0, t = 0, f0 = 0, df0 = 0;
    
    // Ask the user to input the initial guess. Also, ask for the precision of the solution.
    printf("Please input your initial guess x0 and your desired precision: ");
    scanf("%lf %d", &x0, &precision);
    
    printf("\nIteration\tx0\tx1\tf0\tf'0");
    printf("\n___________________________________________________________________\n");
    
    // Newton's method loop.
    do {
        i++;
        t = x1;
        f0 = f(x0);
        df0 = df(x0);
        x1 = x0 - (f0 / df0);     // Newton's method.
        printf("%d %.*lf %.*lf %.*lf %.*lf\n", i, precision + 1, x0, precision + 1, x1, precision + 1, f0, precision + 1, df0, precision + 1);
        x0 = x1;
    } while(fabs(t - x0) > pow(10, -precision));  // Looping until the desired precision is achieved.
    
    printf("___________________________________________________________________\n\n");
    printf("Required root = %.*lf, to %d decimal place(s) precision.\n", precision, t, precision);
    printf("Number of iteration = %d.\n\n", i);
    
    return 0;
}

// The function f(x) = x^3 + 1.5x^2 - 5.75x + 4.37.
double f(double x) {
    return pow(x, 3) + 1.5 * pow(x, 2) - 5.75 * x + 4.37;
}

// The derivative of f(x) --- f'(x).
double df(double x) {
    return 3 * pow(x, 2) + 3 * x - 5.75;
}

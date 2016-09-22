// Q.1. False position method on finding the real roots of x^3 + 1.5x^2 - 5.75x + 4.37 = 0.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double f(double x);

int main () {
    double x0 = 0, x1 = 0, x2 = 0, f0 = 0, f1 = 0, f2 = 0;
    int precision = 0, i = 0;
    
    // Ask the user to input the initial guesses (with f(x0) > 0 and f(x1) < 0). Also, ask for the precision of the solution.
    do {
        printf("Please enter your initial guess x0, such that f(x0) > 0: ");
        scanf("%lf", &x0);
    } while (f(x0) > 0);
    do {
        printf("Please enter your initial guess x1, such that f(x1) < 0: ");
        scanf("%lf", &x1);
    } while (f(x1) < 0);
    
    printf("Please enter your desired precision: ");
    scanf("%d", &precision);
    
    printf("\nIteration\tx0\tx1\tx2\tf0\tf1\tf2");
    printf("\n___________________________________________________________________\n");
    
    do {
        // The loop for false position.
        i++;
        f0 = f(x0);    
        f1 = f(x1);
        x2 = x0 - f0 * (x1 - x0) / (f1 - f0);
        f2 = f(x2);
        
        printf("%d %.*lf %.*lf %.*lf %.*lf %.*lf %.*lf\n", i, precision + 1, x0, precision + 1, x1, precision + 1, x2, precision + 1, f0, precision + 1, f1, precision + 1, f2);
    
        // Selection criterion.
        if (f0 * f2 < 0) {
            x1 = x2;
        }
        else {
            x0 = x2;
        }
    } while(fabs(f2) > pow(10, -precision));  // Looping until the desired precision is achieved.
    
    printf("___________________________________________________________________\n\n");
    printf("Required root = %.*lf, to %d decimal place(s) precision.\n", precision, x2, precision);
    printf("Number of iteration = %d.\n\n", i);
    
    return 0;
}

// The function f(x) = x^3 + 1.5x^2 - 5.75x + 4.37.
double f(double x) {
    return pow(x, 3) + 1.5 * pow(x, 2) - 5.75 * x + 4.37;
}

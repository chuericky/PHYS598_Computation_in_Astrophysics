// Q.2. Perturbation method on finding the small root of x^2 - 10000x + 1 = 0.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double f(double a, double b, double c, double x);

int main () {
    double a = 0, b = 0, c = 0, x0 = 0, x1 = 0, root = 0, t = 0;
    int precision = 0, i = 0;
    
    // Input coeffs of x^2, x then constant term, then the desired precision.
    printf("Enter the coefficients of the quadratic equation: ");
    scanf("%lf %lf %lf", &a, &b, &c);
    
    printf("Enter your initial guess x0: ");
    scanf("%lf", &x0);
    
    printf("Enter your desired precision: ");
    scanf("%d", &precision);
    
    printf("\nIteration\tx0\tx1\tf0\tf1");
    printf("\n____________________________________________________________\n");
    
    do {
        // The loop for perturbation.
        i++;
        x1 = (c + a * pow(x0, 2)) / (-b);
        t = f(a, b, c, x1);
        printf("%d\t%.*lf\t%.*lf\t%.*lf\t%.*lf\n", i, precision + 1, x0, precision + 1, x1, precision + 1, f(a,b,c,x0), precision + 1, f(a,b,c,x1));
        x0 = x1;
    } while(fabs(t) > pow(10, -precision));
    
    printf("___________________________________________________________________\n\n");
    printf("Required root = %.*lf, to %d decimal place(s) precision.\n", precision, x0, precision);
    printf("Number of iteration = %d.\n\n", i);
    
    return 0;
}

// Quadratic equation.
double f(double a, double b, double c, double x) {
    return a * pow(x, 2) + b * x + c;
}
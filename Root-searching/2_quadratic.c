// Q.2. Quadratic formula on finding the small root of x^2 - 10000x + 1 = 0.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main () {
    double a = 0, b = 0, c = 0, root = 0;
    int precision = 0, i = 0;
    
    // Input coeffs of x^2, x then constant term, then the desired precision.
    printf("Enter the coefficients of the quadratic equation: ");
    scanf("%lf %lf %lf", &a, &b, &c);
    
    printf("Enter your desired precision: ");
    scanf("%d", &precision);
    
    // The smaller root.
    root = (-b - pow((pow(b, 2) - 4 * a * c), 0.5)) / (2. * a);
    
    printf("%.*lf\n", precision, root);
    
    return 0;
}

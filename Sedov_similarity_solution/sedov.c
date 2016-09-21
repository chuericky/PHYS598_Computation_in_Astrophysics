// Q.2. Sedov Similarity Solution.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double nu1(double g);
double nu2(double g);
double nu3(double g);
double nu4(double g);
double nu5(double g);
double nu6(double g);
double nu7(double g);
double nu8(double g);
double nu9(double g);
double nu10(double g);
double temp1(double g, double u);
double temp2(double g, double u);
double temp3(double g, double u);
double eqn_up(double g, double u, double eta);
double bisec(double x0, double x1, double g, double eta, double eps);
double find_rhop(double u, double g);
double find_Prat(double u, double g);
double find_Meta(double u, double g);
double simpson(int N, double a, double b, double integ[]);

int main() {
    double gamma = 0, eps = 1e-8, u_up = 0, u_low = 0, a = 0, b = 0, eta = 0.00001, eeta = 0, uprime = 0, rhoprime = 0, Pratio = 0, Pprime = 0, Tratio = 0, M_eta = 0;
    int precision = 14, i = 0, N = 0;

    N = 100000; //(floor)(1 / eta);
    
    // Adjust gamma here.
    gamma = 5./3;
    
    // Initial guess of u for bisection. Lower bound, u_low = (g+1)/2g; Upper bound, u_up = max{5(g+1)/2(3g-1), (g+1)/2}.
    u_low = (gamma + 1) / (2 * gamma);
    a = (5 * (gamma + 1) / (2 * (3 * gamma - 1)));
    b = (gamma + 1) / 2.;
    if (a > b) {
        u_up = b;
    }
    else {
        u_up = a;
    }
    
    printf("\n%s\t%s\t%s\t%s\t%s\t%s\n", "eta", "rho/rho_1", "v/v_1", "P/P_1", "T/T_1", "M(eta)");
    printf("------------------------------------------------------------\n\n");
    
    double *integrands = malloc(N * sizeof(double));
    
    for (i = 1; i < N + 1; i++) {
        eeta = eta * i;
        uprime = bisec(u_low, u_up, gamma, eeta, eps);
        rhoprime = find_rhop(uprime, gamma);
        Pratio = find_Prat(uprime, gamma);
        Pprime = Pratio / pow(eeta, 2);
        Tratio = Pratio / rhoprime; // Assume T/T_1 = (P/P_1)(rho_1/rho)
        M_eta = find_Meta(uprime, gamma);
        integrands[i - 1] = (rhoprime * pow(uprime, 2) + Pprime) * pow(eeta, 4);
        
        printf("%.*lf\t%.*lf\t%.*lf\t%.*lf\t%.*lf\t%.*lf\n", 5, eeta, precision, rhoprime, precision, eeta * uprime, precision, Pratio, precision, Tratio, precision, M_eta);
    }
    
    // eta_0, determined by Simpson's rule.
    printf("\n eta_0 = %.*lf\n", precision, pow((1 / ((32 * M_PI / (25 * (pow(gamma, 2) - 1))) *simpson(N, 0, 1, integrands))), 0.2));

    return 0;
}


// Indices nu1 to nu10.
double nu1(double g) {
    return (13 * pow(g,2) - 7 * g + 12)/((2 * g + 1) * (3 * g - 1));
}

double nu2(double g) {
    return -5 * (g - 1) / (2 * g + 1);
}

double nu3(double g) {
    return 3 / (2 * g + 1);
}

double nu4(double g) {
    return (13 * pow(g,2) - 7 * g + 12)/((2 * g + 1) * (3 * g - 1) * (2 - g));
}

double nu5(double g) {
    return 2 / (g - 2);
}

double nu6(double g) {
    return (13 * pow(g,2) - 7 * g + 12)/(5 * (3 * g - 1) * (2 - g));
}

double nu7(double g) {
    return g / (g - 2);
}

double nu8(double g) {
    return 5 * g / (6 * (2 - g));
}

double nu9(double g) {
    return -(13 * pow(g,2) - 7 * g + 12)/(6 * (2 * g + 1) * (2 - g));
}

double nu10(double g) {
    return -5 * g / (2 * (2 * g + 1));
}

// The three brackets.
double temp1(double g, double u) {
    return (5 * (g + 1) - 2 * u * (3 * g - 1)) / (7 - g);
}

double temp2(double g, double u) {
    return (2 * g * u - (g + 1)) / (g - 1);
}

double temp3(double g, double u) {
    return (g + 1 - 2 * u) / (g - 1);
}

// u'.
double eqn_up(double g, double u, double eta) {
    double tem1 = 0, tem2 = 0;
    
    tem1 = pow(temp1(g, u), nu1(g));
    tem2 = pow(temp2(g, u), nu2(g));
    
    return pow(u, 2) * tem1 * tem2 - pow(eta, -5);
}

// Bisection method to find u'.
double bisec(double x0, double x1, double g, double eta, double eps) {
    double x2 = 0, f0 = 0, f1 = 0, f2 = 0;
    int i = 0;
    const int JMAX=100;
    
    for (i = 0; i < JMAX; i++) {
        x2 = (x0 + x1) / 2.;    // Mid-point.
        f0 = eqn_up(g, x0, eta);
        f1 = eqn_up(g, x1, eta);
        f2 = eqn_up(g, x2, eta);
        
        // Selection criterion.
        if (f0 * f2 < 0) {
            x1 = x2;
        }
        else {
            x0 = x2;
        }
    }
    
    return x2;
}

// Given u', find rho'
double find_rhop(double u, double g) {
    double tem1 = 0, tem2 = 0, tem3 = 0;
    
    tem1 = pow(temp1(g, u), nu4(g));
    tem2 = pow(temp2(g, u), nu3(g));
    tem3 = pow(temp3(g, u), nu5(g));
    
    return tem1 * tem2 * tem3;
}

// Given u', find P/P_1.
double find_Prat(double u, double g) {
    double tem1 = 0, tem2 = 0;
    
    tem1 = pow(temp1(g, u), nu6(g));
    tem2 = pow(temp3(g, u), nu7(g));
    
    return tem1 * tem2 * pow(u, 6./5);
}

// Given u', find M(eta).
double find_Meta(double u, double g) {
    double tem1 = 0, tem2 = 0, tem3 = 0, tau = 0;
    
    tem1 = pow(temp1(g, u), nu9(g));
    tem2 = pow(temp2(g, u), nu10(g));
    tem3 = pow(temp3(g, u), nu8(g));
    
    tau = u * tem1 * tem2 * tem3;
    
    return pow(tau, -6./5);
}

// Simpson's rule.
double simpson(int N, double a, double b, double integ[]) {
    int i = 0;
    double h = 0, sum = 0;
    
    h = (b - a) / N;
    
    if (N == 2) {
        sum += (h / 3.) * (integ[0] + 4 * integ[1] + integ[2]);
    }
    
    else {
        double *f = malloc((N - 1) * sizeof(double));
        
        sum += (h / 3.) * (integ[0] + integ[N - 1]);
        for (i = 0; i < N - 1; i++) {
            if (i % 2 == 0) {
                f[i] = (4/3.) * integ[i + 1];
            }
            else {
                f[i] = (2/3.) * integ[i + 1];
            }
        }
        
        for (i = 0; i < N - 1; i++) {
            sum += h * f[i];
        }
        
        free(f);
    }
    
    return sum;
}

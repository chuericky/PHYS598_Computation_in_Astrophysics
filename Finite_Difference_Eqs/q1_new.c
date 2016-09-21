// Q.1. Rarefaction.

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double q (int i, int j, double **u, double **V);
double t_diff (int i, double **R, double **P, double **V, double gamma, int N_space);


int main() {
    
    // Number of grids in time and spatial coordaintes.
    FILE *filename;
    char name_file[150];
    int N_time = 2000, N_space = 1000;
    int i = 0, j = 0;
    
    // u_ratio is defined as u / u_max. Need to change according to the requirement.
    double u_ratio = -10;
    double gamma = 5./3;
    double temp = 0, dx = 0.1;
    
    // Define double arrays for physical quantities (V, u, R, P, E). Firstly loop over all spatial coordinates, then update the each time steps.
    double **V = malloc (N_time * sizeof (double *));
    double **u = malloc (N_time * sizeof (double *));
    double **R = malloc (N_time * sizeof (double *));
    double **P = malloc (N_time * sizeof (double *));
    double **E = malloc (N_time * sizeof (double *));
    
    for (i = 0; i < N_time; i++) {
        V[i] = malloc (N_space * sizeof (double));
        u[i] = malloc (N_space * sizeof (double));
        R[i] = malloc (N_space * sizeof (double));
        P[i] = malloc (N_space * sizeof (double));
        E[i] = malloc (N_space * sizeof (double));
        for (j = 0; j < N_space; j++) {
            V[i][j] = 0;
            u[i][j] = 0;
            R[i][j] = 0;
            P[i][j] = 0;
            E[i][j] = 0;
        }
    }
    
    // Initial values for all physical quantities. All quantities are normalized.
    for (j = 1; j < N_space; j++) {
        // Values at t = 0.
        V[0][j] = 1;
        u[0][j] = 0;
        P[0][j] = 1;
    }
    
    for (i = 0; i < N_time - 1; i++) {
        // Values at the piston boundary.
        V[i][0] = 1;
        u[i][0] = ((2 * sqrt(gamma)) / (gamma - 1)) * u_ratio;
        P[i][0] = 1;
        u[i][N_space - 1] = 0;
        V[i][N_space - 1] = 1;
        P[i][N_space - 1] = 1;
    }
    
    // At time t = 0:
    // Spatial coordinates and energy at time = 0 at each spatial coordinates.
    for (j = 0; j < N_space; j++) {
        R[0][j] = j * dx;
        E[0][j] = (P[0][j] * V[0][j]) / (gamma - 1);    // Equation of State.
    }
    
    
    // Update every time step and spatial coordinates. For t > 0 and x > 0.
    for (i = 0; i < N_time - 1; i++) {  // N_time - 1
        if (i % 100 == 0) printf("%d\n", i);
        
        // Step 1: u
        for (j = 1; j < N_space; j++) {
            temp = t_diff(i, R, P, V, gamma, N_space) * (2 * (P[i][j] + q(i, j, u, V) - q(i, j - 1, u, V) - P[i][j - 1]) / (R[0][j + 1] - R[0][j - 1])) * (V[0][j] + V[0][j - 1]) / 2.;
            u[i + 1][j] = u[i][j] - temp;
                   
                   
                   //c = sqrt(gamma * P[i][j] * V[i][j]);
                   //temp = (R[i][j + 1] - R[i][j]) / c;)
            //printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\n\n", t_diff(i, R, P, V, gamma, N_space), R[0][j], P[i][j], u[i][j], V[i][j], q(i, j, u, V), temp);
        }
        
        // Step 2: R
        for (j = 0; j < N_space; j++) {
            R[i + 1][j] = R[i][j] + t_diff(i, R, P, V, gamma, N_space) * u[i + 1][j];
            //printf("%f\t%f\t%f\t%f\n", R[i + 1][j], R[i][j], t_diff(i, R, P, V, gamma, N_space), u[i + 1][j]);
        }
        
        // Step 3: V
        for (j = 0; j < N_space; j++) {
            V[i + 1][j] = V[0][j] * (R[i + 1][j + 1] - R[i + 1][j]) / (R[0][j + 1] - R[0][j]);
            //printf("%f\t%f\t%f\t%f\t%f\t%f\n", V[i + 1][j], V[0][j], R[i + 1][j + 1], R[i + 1][j], R[0][j + 1], R[0][j]);
        }
        
        // Step 4 and 5: E, P
        for (j = 0; j < N_space; j++) {
            E[i + 1][j] = (E[i][j] - ((P[i][j] / 2) + q(i + 1, j, u, V)) * (V[i + 1][j] - V[i][j])) / (1 + 0.5 * (gamma - 1) * (V[i + 1][j] - V[i][j]) / V[i + 1][j]);
            P[i + 1][j] = (gamma - 1) * E[i + 1][j] / V[i + 1][j];
            //printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\n", E[i + 1][j], E[i][j], P[i][j], q(i + 1, j, u, V), V[i + 1][j], V[i][j], P[i + 1][j]);
        }
    }
    
    for (i = 500; i < 501; i++) { // i = 0; N_time - 1
        //printf("j\tR\t1/V\tE\tP\tu\tdt\tq\n");
        strcpy (name_file, "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps5/data/10_500.dat");
        
        filename = fopen (name_file, "w");
        
        for (j = 0; j < N_space; j++) {
            //printf("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n", j, R[i][j], 1/V[i][j], E[i][j], P[i][j], u[i][j], t_diff(i, R, P, V, gamma, N_space), q(i, j, u, V));
            fprintf(filename, "%d\t%f\t%f\t%f\t%f\n", j, R[i][j], 1/V[i][j], P[i][j], u[i][j]);
            
        }
        
        fclose (filename);
        
        printf("j\tR\t1/V\tE\tP\tu\tdt\tq\n");
    }
    
    
    
    
    
    free (V);
    free (u);
    free (R);
    free (P);
    free (E);
    
    
    
    return 0;
}



double q (int i, int j, double **u, double **V) {
    double a = 1.5;
    
    // Make sure i >= 0, and j >= 0.
    if (i < 0) {
        printf("Error! i < 0.\n");
        exit(1);
    }
    if (j < 0) {
        printf("Error! j < 0.\n");
        exit(1);
    }
    
    if (u[i][j + 1] - u[i][j] < 0) {
        if (i == 0) {
            return a * a * pow((u[i][j + 1] - u[i][j]), 2) / V[i][j];
        }
        else {
            return 2 * a * a * pow((u[i][j + 1] - u[i][j]), 2) / (V[i][j] + V[i + 1][j]);
        }
    }
    else    return 0;
}


double t_diff (int i, double **R, double **P, double **V, double gamma, int N_space) {
    // Time step.
    double b = 0.2, min = 9999999999., c = 0, temp = 0;
    int j = 0;
    
    for (j = 0; j < N_space - 1; j++) {
        c = sqrt(gamma * P[i][j] * V[i][j]);
        temp = (R[i][j + 1] - R[i][j]) / c;
        //printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\n", P[i][j], V[i][j], c, R[i][j + 1], R[i][j], temp, min);
        if (temp < min) min = temp;
    }
    if (min < 0) {
        printf("Error! t_diff < 0.\n");
        exit(1);
    }
    return b * min;
}
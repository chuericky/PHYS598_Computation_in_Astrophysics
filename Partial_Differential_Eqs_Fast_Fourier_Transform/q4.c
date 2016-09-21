// Q.4. Fourier Transform of Gaussian function h(t) = exp(-t^2).

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "math.h"

int main() {
    FILE *foufile, *invfoufile;
    char name_fou[150], name_invfou[150], steps[10];
    double t0 = -10.0, t1 = 10.0, t = 0, f = 0;
    int N = pow(2, 5);        // N, number of samples.
    float *data;
    int i = 0;
    
    strcpy (name_fou, "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps10/q_4/data/");
    sprintf (steps, "%d", N);
    strcat (name_fou, steps);
    strcat (name_fou, "_fourier.dat");
    
    foufile = fopen (name_fou, "w");
    
    strcpy (name_invfou, "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps10/q_4/data/");
    sprintf (steps, "%d", N);
    strcat (name_invfou, steps);
    strcat (name_invfou, "_invfourier.dat");
    
    invfoufile = fopen (name_invfou, "w");
    
    data = vector(1, 2 * N + 1);
    
    for (i = 1; i <= 2 * N - 1; i++) {
        t = t0 + i * (t1 - t0) / (2 * N);
        data[i] = exp(-pow(t, 2));
        data[i + 1] = 0;
        //printf("%d\t%E\t%E\n", i, data[i], data[i + 1]);
    }
    
    // Fourier Transform
    four1(data, N, 1);
    
    for (i = N/2 + 1; i <= N; i++) {
        f = (i - 1 - N) / (t1 - t0);
        fprintf(foufile, "%d\t%f\t%E\t%E\n", (i - 1), f, data[2*i - 1] * (t1 - t0) / N, data[2*i] * (t1 - t0) / N);
    }
    for (i = 1; i < N/2 + 1; i++) {
        f = (i - 1) / (t1 - t0);
        fprintf(foufile, "%d\t%f\t%E\t%E\n", (i - 1), f, data[2*i - 1] * (t1 - t0) / N, data[2*i] * (t1 - t0) / N);
    }
    
    // Inverse Fourier Transform
    four1(data, N, -1);
    
    for (i = 1; i <= N; i++) {
        t = t0 + i * (t1 - t0) / N;
        fprintf(invfoufile, "%d\t%f\t%E\t%E\n", (i - 1), t, data[2*i - 1] / N, data[2*i] / N);  // Have to divide N manually.
    }
    
    
    free_vector(data, 1, 2 * N + 1);
    
    fclose (foufile);
    fclose (invfoufile);
    
    return 0;
}
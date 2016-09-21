// Q.6. Pendulum.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "math.h"

int main() {
    FILE *foufile_x, *foufile_v, *penfile;
    char name_fou_x[150], name_fou_v[150], name_pen[150], steps[10];
    int N = pow(2, 14);        // N, number of samples.
    float *t, *x, *v;
    float dt = 0.1, power_x = 0, power_v = 0, f = 0;
    int i = 0, count = 0, ps6_N = 1001;
    
    strcpy (name_pen, "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps10/q_6/data/2.dat");
    
    penfile = fopen (name_pen, "r");
    
    strcpy (name_fou_v, "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps10/q_6/data/2_pen_v.dat");
    
    foufile_v = fopen (name_fou_v, "w");
    
    strcpy (name_fou_x, "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps10/q_6/data/2_pen_x.dat");
    
    foufile_x = fopen (name_fou_x, "w");
    
    t = vector(1, 2 * N + 1);
    x = vector(1, 2 * N + 1);
    v = vector(1, 2 * N + 1);
    
    float *t_arr = malloc (ps6_N * sizeof (float));
    float *x_arr = malloc (ps6_N * sizeof (float));
    float *v_arr = malloc (ps6_N * sizeof (float));
    
    count = 0;
    
    while (fscanf(penfile, "%f %f %f", &t_arr[count], &x_arr[count], &v_arr[count])== 3)
        count++;
    
    for (i = 0; i < 2 * N + 1; i++) {
        t[i] = 0;
        x[i] = 0;
        v[i] = 0;
    }
    
    for (i = 1; i < ps6_N + 1; i++) {
        x[2 * i - 1] = x_arr[i - 1];
        v[2 * i - 1] = v_arr[i - 1];
    }
    
    // Fourier Transform
    four1(v, N, 1);
    four1(x, N, 1);
    
    for (i = 1; i < N/2 + 1; i++) {
        f = (i - 1) / (N * dt);
        power_v = 2 * (pow(v[2*i - 1], 2) + pow(v[2*i], 2));
        fprintf(foufile_v, "%d\t%f\t%f\t%E\t%E\t%E\n", (i - 1), f, 2 * M_PI * f / (2./3), v[2*i - 1], log(power_v), power_v);
        power_x = 2 * (pow(x[2*i - 1], 2) + pow(x[2*i], 2));
        fprintf(foufile_x, "%d\t%f\t%f\t%E\t%E\t%E\n", (i - 1), f, 2 * M_PI * f / (2./3), x[2*i - 1], log(power_x), power_x);
    }
    
    free_vector(t, 1, 2 * N + 1);
    free_vector(x, 1, 2 * N + 1);
    free_vector(v, 1, 2 * N + 1);
    
    free(t_arr);
    free(x_arr);
    free(v_arr);
    
    fclose (foufile_x);
    fclose (foufile_v);
    fclose (penfile);
    
    return 0;
}
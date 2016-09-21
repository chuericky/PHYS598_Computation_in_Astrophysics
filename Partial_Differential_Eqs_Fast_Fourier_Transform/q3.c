// Q.3. Linear Elliptic PDE by multigrid method.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "nr.h"
#include "nrutil.h"
#include "math.h"

int main() {
    FILE *sorufile, *sorFfile;
    char name_sor_u[150], name_sor_F[150];
    double **f, **u, t;
    int i = 0, j = 0, midl = 0;
    const int Nstep = 4, JMAX = 33;
    time_t start = 0, end = 0;
    
    strcpy (name_sor_u, "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps10/q_3/data/u.dat");
    strcpy (name_sor_F, "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps10/q_3/data/F_float.dat");
    sorufile = fopen (name_sor_u, "w");
    sorFfile = fopen (name_sor_F, "w");
    
    midl = (JMAX + 1) / 2;
    
    f = dmatrix(1, JMAX + 1, 1, JMAX + 1);
    u = dmatrix(1, JMAX + 1, 1, JMAX + 1);
    
    for (i = 1; i <= JMAX; i++) {
        for (j = 1; j <= JMAX; j++) {
            f[i][j] = 0;
            u[i][j] = 0;
        }
    }
    u[midl][midl] = 2.;    // Setting the density at everywhere.
    
    // SOR method.
    start = clock();
    mglin(u, JMAX, 10);
    end = clock();
    t=(end - start)/CLOCKS_PER_SEC;
    
    printf("Time taken %E seconds.\n", t);
                           
    for (i = 1; i <= JMAX; i+=Nstep) {
        for (j = 1; j <= JMAX; j+=Nstep) {
            fprintf(sorufile, "%.*f\t", 5, u[i][j]);
        }
        fprintf(sorufile, "\n");
    }
    for (i = Nstep + 1; i < JMAX; i+=Nstep) {
        for (j = Nstep + 1; j < JMAX; j+=Nstep) {
            f[i][j] = u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] - 4 * u[i][j];
            fprintf(sorFfile, "%.*f\t", 5, f[i][j] * pow((JMAX - 1), 2));
        }
        fprintf(sorFfile, "\n");
    }
    
    free_dmatrix(f, 1, JMAX + 1, 1, JMAX + 1);
    free_dmatrix(u, 1, JMAX + 1, 1, JMAX + 1);
    
    fclose (sorufile);
    fclose (sorFfile);
    
    return 0;
}
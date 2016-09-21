// Q.2. Non-Linear Elliptic PDE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "nr.h"
#include "nrutil.h"
#include "math.h"

int main() {
    FILE *sorufile, *sorFfile;
    char name_sor_u[150], name_sor_F[150];
    double **a, **b, **c, **d, **e, **f, **u;
    int i = 0, j = 0, midl = 0;
    const int Nstep = 4, JMAX = 33;
    double rjac = cos(M_PI / JMAX);
    
    strcpy (name_sor_u, "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps10/q_2/data/u.dat");
    strcpy (name_sor_F, "/Users/rickyccy/Documents/Urbana-Champaign/Courses/PHYS598/ps10/q_2/data/F_float.dat");
    sorufile = fopen (name_sor_u, "w");
    sorFfile = fopen (name_sor_F, "w");
    
    midl = (JMAX + 1) / 2;
    
    a = dmatrix(1, JMAX + 1, 1, JMAX + 1);
    b = dmatrix(1, JMAX + 1, 1, JMAX + 1);
    c = dmatrix(1, JMAX + 1, 1, JMAX + 1);
    d = dmatrix(1, JMAX + 1, 1, JMAX + 1);
    e = dmatrix(1, JMAX + 1, 1, JMAX + 1);
    f = dmatrix(1, JMAX + 1, 1, JMAX + 1);
    u = dmatrix(1, JMAX + 1, 1, JMAX + 1);
    
    for (i = 1; i <= JMAX; i++) {
        for (j = 1; j <= JMAX; j++) {
            a[i][j] = 1;
            b[i][j] = 1;
            c[i][j] = 1;
            d[i][j] = 1;
            e[i][j] = -4;
            f[i][j] = 0;
            u[i][j] = 0;
        }
    }
    f[midl][midl] = 2./(pow((JMAX - 1), 2));    // Setting the density at everywhere.
    
    // SOR method.
    sor(a, b, c, d, e, f, u, JMAX, rjac);
                           
    for (i = 1; i <= JMAX; i+=Nstep) {
        for (j = 1; j <= JMAX; j+=Nstep) {
            fprintf(sorufile, "%.*f\t", 5, u[i][j]);
        }
        fprintf(sorufile, "\n");
    }
    for (i = Nstep + 1; i < JMAX; i+=Nstep) {
        for (j = Nstep + 1; j < JMAX; j+=Nstep) {
            f[i][j] = u[i + 1][j] + u[i - 1][j] + u[i][j + 1] + u[i][j - 1] - 4 * u[i][j] + pow(u[i][j] / JMAX, 2);
            fprintf(sorFfile, "%.*f\t", 5, f[i][j]);
        }
        fprintf(sorFfile, "\n");
    }
    
    free_dmatrix(a, 1, JMAX + 1, 1, JMAX + 1);
    free_dmatrix(b, 1, JMAX + 1, 1, JMAX + 1);
    free_dmatrix(c, 1, JMAX + 1, 1, JMAX + 1);
    free_dmatrix(d, 1, JMAX + 1, 1, JMAX + 1);
    free_dmatrix(e, 1, JMAX + 1, 1, JMAX + 1);
    free_dmatrix(f, 1, JMAX + 1, 1, JMAX + 1);
    free_dmatrix(u, 1, JMAX + 1, 1, JMAX + 1);
    
    fclose (sorufile);
    fclose (sorFfile);
    
    return 0;
}
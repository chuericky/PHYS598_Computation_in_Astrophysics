#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main() {

    FILE *filename1,*filename2,*filename3;
    char file1[150],file2[150],file3[150];
    int N_time = 1001;
    int N_space = 1001;
    int i = 0, j = 0;
    
    float x0 = -10.0, x1 = 10.0;
    double *r, *u;
    
    double dx = (x1 - x0) / (N_space - 1.);     // Step in x.
    double dt = 0.1 * dx * dx;                  // Time step.
    double t= 0.0;

    u = new double[N_space];
    r = new double[N_space]; //values of u from previous time step

    // Initial condition.
    for (i = 0; i < N_space; i++) {
        u[i] = exp(-pow((x0 + i * dx), 2));
    }
       
    // Finite difference method    
    for (j = 1; j < N_time; j++) {	

	for(i=0; i < N_space; i++){
		if (i == 0 || i == N_space - 1) r[i] = 0; //Boundary Conditions
         	else    r[i] = u[i];//store the previous u information
        }					
		
	// Solve for other u's
	for(i=1; i<N_space-1; i++){
		u[i]=r[i] + dt/(dx*dx) * (r[i+1]-2*r[i]+r[i-1]);
	}
				
	t+=dt; 

	// Print out the information at a specified timestep                                
        if (j == 100) {
	    strcpy (file1, "/home/quantum-monkey/workspace/CPAcodes/ps9/data/p1data1.dat");    
            filename1 = fopen (file1, "w"); 
            for (i = 0; i < N_space; i++) {
                fprintf(filename1,"%d\t%d\t%f\t%f\n", j, i, x0 + i * dx, u[i]);
            }
	fclose (filename1);
        }

	if (j == 700) {
	    strcpy (file2, "/home/quantum-monkey/workspace/CPAcodes/ps9/data/p1data2.dat");    
            filename2 = fopen (file2, "w"); 
            for (i = 0; i < N_space; i++) {
                fprintf(filename2,"%d\t%d\t%f\t%f\n", j, i, x0 + i * dx, u[i]);
            }
	fclose (filename2);
        }

	if (j == 900) {
	    strcpy (file3, "/home/quantum-monkey/workspace/CPAcodes/ps9/data/p1data3.dat");    
            filename3 = fopen (file3, "w"); 
            for (i = 0; i < N_space; i++) {
                fprintf(filename3,"%d\t%d\t%f\t%f\n", j, i, x0 + i * dx, u[i]);
            }
	fclose (filename3);
        }
    }    
         
    free(u);
    free(r);
    
    return 0;
}


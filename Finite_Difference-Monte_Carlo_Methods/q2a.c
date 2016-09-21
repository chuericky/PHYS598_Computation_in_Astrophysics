#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


int main() {

    FILE *filename1,*filename2,*filename3,*filename4;
    char file1[150],file2[150],file3[150],file4[150];
    int N_time = 10001;
    int N_space = 1001;
    int i = 0, j = 0;
    
    //float x0 = -10.0, x1 = 10.0;
    float x0 = 0.0, x1 = 100.0;
    double *r, *u, *E,*E0 ;
    
    double dx = (x1 - x0) / (N_space - 1.);     // Step in x.
    double dt = 0.1 * dx * dx;                  // Time step.
    double t= 0.0;

    u = new double[N_space];
    r = new double[N_space]; //values of u from previous time step
    E = new double[N_time]; 
    E0 = new double[N_time];
    
    // Initial condition.
    for (i = 0; i < N_space; i++) {
	u[i] = 1.0;
    }
       
    // Finite difference method    
    for (j = 1; j < N_time; j++) {	        

	for(i=0; i < N_space; i++){
		r[i] = u[i];//store the previous u information
        }
		
	// Solve for u[1] in new time step
	u[1]=r[1] + dt/(dx*dx) * (r[2]-2*r[1]+r[0]);
		
	// Go back for the new u[0] (Based on b.c.)
	u[0]=u[1]/(1+dx);
				
	// Solve for other u's
	for(i=2; i<N_space-1; i++){
		u[i]=r[i] + dt/(dx*dx) * (r[i+1]-2*r[i]+r[i-1]);
	}
				
	t+=dt; 

	// Energy conservation
	for (i=0; i < N_space; i++){
		E[j] += u[i]*dx ;
	}
	E0[j] = u[0]; //RHS of integral
	

	// Print out the information at a specified timestep                                
        if (j == 100) {
	    strcpy (file1, "/home/quantum-monkey/workspace/CPAcodes/ps9/data/p2data1.dat");    
            filename1 = fopen (file1, "w"); 
            for (i = 0; i < N_space; i++) {
                fprintf(filename1,"%d\t%d\t%f\t%f\n", j, i, x0 + i * dx, u[i]);
            }
	fclose (filename1);
        }

	if (j == 5000) {
	    strcpy (file2, "/home/quantum-monkey/workspace/CPAcodes/ps9/data/p2data2.dat");    
            filename2 = fopen (file2, "w"); 
            for (i = 0; i < N_space; i++) {
                fprintf(filename2,"%d\t%d\t%f\t%f\n", j, i, x0 + i * dx, u[i]);
            }
	fclose (filename2);
        }

	if (j == 7000) {
	    strcpy (file3, "/home/quantum-monkey/workspace/CPAcodes/ps9/data/p2data3.dat");    
            filename3 = fopen (file3, "w"); 
            for (i = 0; i < N_space; i++) {
                fprintf(filename3,"%d\t%d\t%f\t%f\n", j, i, x0 + i * dx, u[i]);
            }
	fclose (filename3);
        }
    }

    //Print out energy conservation
    strcpy (file4, "/home/quantum-monkey/workspace/CPAcodes/ps9/data/p2energy.dat");    
    filename4 = fopen (file4, "w"); 
    for ( j = 1; j < N_time; j++) {
	if (j%50 ==0) 
		fprintf(filename4,"%d\t%f\t%f\n",j, j*dt, (E[j]-E[j-1]+E0[j]*dt)/(E[j]-E[j-1]));
    }

    fclose(filename4);

    free(u);
    free(r);
    free(E);
    free(E0);

    return 0;
}


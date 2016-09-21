#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double t_diff (int j,double *u, int N_space, double dx);


int main() {

    FILE *filename1,*filename2,*filename3;
    char file1[150],file2[150],file3[150];
    int N_time = 50001;
    int N_space = 101;
    int i = 0, j = 0;
    
    float x0 = -1.0, x1 = 1.0;
    double *r, *u;
    
    double dx = (x1 - x0) / (N_space - 1.);     // Step in x.
    double *dt;                  // Time step
    double t= 0.0;
    dt = new double[N_time];

    u = new double[N_space];
    r = new double[N_space]; //values of u from previous time step
    
    // Initial condition.
    for (i = 0; i < N_space; i++) {
	u[i]= x0 + i * dx ;
    }
       
    // Finite difference method    
    for (j = 1; j < N_time; j++) {	        
	dt[j] = t_diff(j,u,N_space,dx); //determine time step
	//printf("%.*f",10,t_diff(j,u,N_space,dx));
	for(i=0; i < N_space; i++){
	    if (i==0) r[i] = -1.0;
	    else if (i == N_space -1) r[i] = 1.0;	//Boundary condition	
	    else r[i] = u[i];
        }
						
	// Solve for other u's
	for(i=1; i<N_space-1; i++){
		u[i]=r[i] + dt[j]/(dx*dx) * (pow(r[i+1],3.0)-2*pow(r[i],3.0)+pow(r[i-1],3.0));
	}
	
	//u[0] = -1.0;
	//u[N_space] = 1.0;
 			
	t+=dt[j]; 
	//printf("%.*f",10,t);	

	// Print out the information at a specified timestep                                
        if (j == 100) {
	    strcpy (file1, "/home/quantum-monkey/workspace/CPAcodes/ps9/data/p3data1.dat");    
            filename1 = fopen (file1, "w"); 
            for (i = 0; i < N_space; i++) {
                fprintf(filename1,"%d\t%d\t%f\t%f\n", j, i, x0 + i * dx, u[i]);
            }
	fclose (filename1);
        }

	if (j == 5000) {
	    strcpy (file2, "/home/quantum-monkey/workspace/CPAcodes/ps9/data/p3data2.dat");    
            filename2 = fopen (file2, "w"); 
            for (i = 0; i < N_space; i++) {
                fprintf(filename2,"%d\t%d\t%f\t%f\n", j, i, x0 + i * dx, u[i]);
            }
	fclose (filename2);
        }

	if (j == 50000) {
	    strcpy (file3, "/home/quantum-monkey/workspace/CPAcodes/ps9/data/p3data3.dat");    
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

    // Time step selection ( to ensure stability) 
double t_diff (int j,double *u, int N_space, double dx) {

    double b = 0.1, min = 9999999., temp = 0;
    int k = 0;
    
    for (k = 0; k < N_space; k++) {
        temp = pow(u[k], -2.0)*dx*dx/3.0;	
        if (temp < min) min = temp;
    }    
    if (min < 0) {
        printf("Error! t_diff < 0.\n");
        exit(1);
    }
    //printf("%d",b*min);
    return b * min;
}



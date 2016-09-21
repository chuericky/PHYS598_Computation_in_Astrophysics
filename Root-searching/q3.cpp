// Q.3. 

#include "nr.h"
#include <iostream>
#include <math.h>
#include <iomanip>
#include <vector>

int main()
{	
	int ntrial = 20;
	Vec_DP xtry(2);
	double tolx = 1e-8, tolf = 1e-8;
	
	xtry[0] = 0.1;
	xtry[1] = 1.5;
	
	NR::mnewt(ntrial, xtry, tolx, tolf);
	
	cout << setprecision(7) << "The solution is (x,y) = (" << xtry[0] << "," << xtry[1] << ")." << endl;
	
	return 0;
}

void usrfun(Vec_I_DP &x, Vec_O_DP &fvec, Mat_O_DP &fjac)
{
	fvec[0] = x[0]*x[0]+x[1]*x[1]-1;
	fvec[1] = x[0]-x[1];
	fjac[0][0] = 2*x[0];
	fjac[0][1] = 2*x[1];
	fjac[1][0] = 1;
	fjac[1][1] = -1;
}
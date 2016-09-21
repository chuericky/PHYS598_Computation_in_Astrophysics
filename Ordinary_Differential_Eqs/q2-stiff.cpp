// Q.2. Stiff ODE integrator

#include <iostream>
#include <math.h>
#include <iomanip>
#include <cmath>
#include "nr.h"
using namespace std;

void derivs(const DP x, Vec_I_DP &y, Vec_O_DP &dydx);

int kmax, kount;        // defining declarations
DP dxsav;
Vec_DP *xp_p;
Mat_DP *yp_p;

int main() {
	const int NP = 3;
    DP eps, hstart, x1 = 0., x2 = 1000.;
    Vec_DP y(NP);
    int nbad, nok;
    
    eps = 0.005;
    hstart = 0.005;
    kmax = 0;
          
    // initial conditions
    y[0] = 1.0;
    y[1] = 1.0;
    y[2] = 0.0;
          
    NR::odeint(y, x1, x2, eps, hstart, 0.0, nok, nbad, derivs, NR::stiff);
    
	return 0;	
}


void derivs(const DP x, Vec_I_DP &y, Vec_O_DP &dydx){
    dydx[0] = -0.013 * y[0] - 1000 * y[0] * y[2];
    dydx[1] = -2500 * y[1] * y[2];
    dydx[2] = -0.013 * y[0] - 1000 * y[0] * y[2] - 2500 * y[1] * y[2];
    
}
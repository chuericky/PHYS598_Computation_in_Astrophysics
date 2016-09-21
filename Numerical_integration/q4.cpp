#include <iostream>
#include <math.h>
#include <iomanip>

using namespace std;

double findy(double x);
double simpson(double a, double b, int Ndiv, double epsilon);

int main()
{
	double delta;
	double a, b;
	int Ndiv=10;
	double integral, integral_old=-99999, epsilon=1e-8;
	
	delta=0.1;
	a=delta;
	b=1-delta;
	
	while(fabs((integral-integral_old)/integral_old)>epsilon)
	{
		integral = 3*pow(delta,1.0/3.0)+1.0/4.0*pow(delta,4.0/3.0) + simpson(a,b,Ndiv,epsilon) + 3.0/2.0*pow(delta,2.0/3.0)+2.0/5.0*pow(delta,5.0/3.0);
		
		cout << setprecision(9) << "delta = " << delta << "\t integral = " << integral << endl;
		delta/=2;
		a=delta;
		b=1-delta;
	}	
	return 0;
}

double findy(double x)
{
	return pow(x,-2.0/3.0)*pow(1-x,-1.0/3.0);
}

double simpson(double a, double b, int Ndiv, double epsilon)
{
	double pi=4*atan(1.0);
	double h, x, sum, sumold=-99999, hnew;
	int i, count=1;	
	
	//cout << "Simpsons Rule" << endl;
	
	h=(b-a)/Ndiv;
	
	sum=findy(a)+findy(b);
	x=a+h;
	while(x<b)
	{
		sum+=2*findy(x);
		x+=h;	
	}
	x=a+h/2;
	while(x<b)
	{
		sum+=4*findy(x);
		x+=h;	
	}
	
	sum*=(h/2)/3;
	
	//cout << setprecision(11) << "iteration = " << count << "\t integral = " << sum << endl;
	count++;
	
	while(fabs((sum-sumold)/sumold)>epsilon)
	{
		sumold=sum;
		sum/=(h/2)/3;
		h/=2;
		x=a+h;
		while(x<b)
		{
			sum-=2*findy(x);
			x+=2*h;
		}
		x=a+h/2;
		while(x<b)
		{
			sum+=4*findy(x);
			x+=h;
		}
		sum*=(h/2)/3;
		
		//cout << setprecision(11) << "iteration = " << count << "\t integral = " << sum << endl;
		count++;
	}
	return sum;	
}
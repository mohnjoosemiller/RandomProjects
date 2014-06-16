/*
	Heath Henley 
	University of Rhode Island

	Simpson's Rule and simple Monte Carlo integration
*/
#include <iostream> 
#include <vector>
#include <stdlib.h> 
using namespace std; 


// function to integrate f(x) = x^2
double f(double x ) { return x*x; }
// function to integrate g(x) = x^10
double g(double x ) { return pow(x,10); }
// function to integrate f(X) = x1^4*...x10^4 * exp( x1+ ... + x10 ) 
double f10D(vector<double> &x ) 
{ 
	double f = 1; 
	for (int i = 0; i < 10; i++) 
		f*= pow(x[i], 4) * exp( x[i] );
	return f; 
}
double rho(double x){return exp(-x*x);}

double rho(vector<double> &x)
{
	double rho = 1; 
	for ( int i = 0; i < 10; i++) 
		rho*= exp(-x[i]*x[i]);
	return rho;
}

double f(vector<double> &x)
{
	double f = 1; 
	for ( int i = 0; i < 10; i++) 
		f*= x[i]*x[i];
	return f;
}

void simpsonsRule(double a , double b , int points)
{
	
	// step size 
	double h = abs(b - a)/points;

	a  = 0; b = 1; 

	double I = 0; 
	// interior points
	for ( int i = 1; i < points-1; i++) 
	{	
		I += f(a+i*h) ;
	} 
	// end points 
	I += f(a)/2.0; 
	I += f(b)/2.0;
	I *= h; 

	cout << "The value of the integral using simpson's rule is " << I << endl;
}
void simpsonsRule2(double a , double b , int points)
{
	
	// step size 
	double h = abs(b - a)/points;

	a  = 0; b = 1; 

	double I = 0; 
	// interior points
	for ( int i = 1; i < points-1; i++) 
	{	
		I += g(a+i*h) ;
	} 
	// end points 
	I += g(a)/2.0; 
	I += g(b)/2.0;
	I *= h; 

	cout << "The value of the integral using simpson's rule is " << I << endl;
}
void monteCarlo(int points)
{
	double I = 0,I2 = 0; 
	
	for ( int i = 0; i < points; i++)
	{
		double x = (double)rand()/RAND_MAX;
		I  += f(x);
		I2 += f(x)*f(x);
	}
	I  *= (1.0/points);
	I2 *= (1.0/points); 
	double std = sqrt(I2-I*I)/(points-1);
	cout << "The value of the integral using Monte Carlo is " << I << " +- " << 2*std<<  endl;
}
void monteCarlo2(int points)
{
	double I = 0,I2 = 0; 
	
	for ( int i = 0; i < points; i++)
	{
		double x = (double)rand()/RAND_MAX;
		I  += g(x);
		I2 += g(x)*g(x);
	}
	I  *= (1.0/points);
	I2 *= (1.0/points); 
	double std = sqrt(I2-I*I)/(points-1);
	cout << "The value of the integral using Monte Carlo is " << I << " +- " << 2*std<<  endl;
}
void monteCarlo10D(int points)
{
	double I = 0,I2 = 0; 
	vector<double> x(10,0);
	for ( int i = 0; i < points; i++)
	{
		// assign random values to each x1 ... x10 
		for ( int j = 0 ; j < 10; j++)
			x[j] = (double) rand()/RAND_MAX;
		I  += f10D(x);
		I2 += f10D(x)*f10D(x);
	}
	I  *= (1.0/points);
	I2 *= (1.0/points); 
	double std = sqrt(I2-I*I)/(points-1);
	cout << "The value of the integral using Monte Carlo is " << I << " +- " << 2*std<<  endl;
}

void monteCarlo_inf(int equil, int prod)
{
	double I = 0,I2 = 0;
	double x = 0, x_new;
	double delta = 50.0, r;
	int accepted = 0; 
	
	for ( int i = 0; i < equil+prod; i++)
	{
		x_new = x + ((double)rand()/RAND_MAX - 0.5)*delta;
		if (rho(x_new)/rho(x) >= ((double)rand()/RAND_MAX))
		{
			x = x_new;
			accepted++;
		}
		I  += f(x);
		I2 += f(x)*f(x);

		if ( i % 100 == 0 && i < equil ) 
		{
			// calculate acceptance ratio
			r = (double)accepted / i;
			if ( r > 0.5 ) 
				delta = delta*1.05;
			else 
				delta = delta*0.95;
		}
		if ( i == equil) 
		{
			cout << "\tSwitch to prod ...." << endl; 
			I = 0.0; 
			I2 = 0.0;
		}
	}
	I  *= (1.0/prod);
	I2 *= (1.0/prod); 
	double std = sqrt(I2-I*I)/(prod-1);
	cout << "The value of the integral using Monte Carlo is " << I << " +- " << 2*std<<  endl;
	cout << "Acceptance ratio =  " << r << ", Delta = " << delta << endl;
	
}

void monteCarlo_inf_10D(int equil, int prod)
{
	double I = 0,I2 = 0;
	double delta = 50, r;
	int accepted = 0; int j = 0; 
	vector < double > x(10,0.0), x_new(10, 0.0);
	
	for ( int i = 0; i < equil+prod; i++)
	{
		j = rand() % 10;
		x_new[j] = x[j] + ((double)rand()/RAND_MAX - 0.5)*delta;

		if (rho(x_new)/rho(x) >= ((double)rand()/RAND_MAX))
		{
			x = x_new;
			accepted++;
		}
		else 
			x_new = x;

		I  += f(x);
		I2 += f(x)*f(x);

		if ( i % 100 == 0 && i < equil) 
		{
			// calculate acceptance ratio
			r = (double)accepted / i;
			if ( r > 0.5 ) 
				delta = delta*1.05;
			else 
				delta = delta*0.95;
		}
		if ( i == equil) 
		{
			cout << "\tSwitch to prod ...." << endl; 
			I = 0.0; 
			I2 = 0.0;
		}
	}
	I  *= (1.0/prod);
	I2 *= (1.0/prod); 
	double std = sqrt(I2-I*I)/(prod-1);
	cout << "The value of the integral using Monte Carlo is " << I << " +- " << 2*std<<  endl;
	cout << "Acceptance ratio =  " << r << ", Delta = " << delta << endl;
}

int main ()
{
	
	double a = 0, b = 1;
	int points = 1e5; 
	int mc_points = points; 
	int equil = 1e6;
	int prod = 2e7;

	/* 
	Finite integration limits 
	*/

	// Use Simpson's rule for f = f(x) = x^2  
	simpsonsRule(a, b, points);

	// Use Simpson's rule for g = g(x) = x^10  
	simpsonsRule2(a, b, points);

	// Monte carlo integrate f(x) = x^2 
	monteCarlo(mc_points);

	// Monte carlo integrate g(x) = x^10 
	monteCarlo2(mc_points);

	// Monte Carlo Integrate f(X) = X^4 * exp (X) - 10 dimensions
	//monteCarlo10D(mc_points);


	/* 
	Integrals over all space 
	*/
	// Monte Carlo Integrate:  Integral( x^2 *exp(-x^2) )/ Integral( exp(-x^2) ) both from -inf to inf 
	monteCarlo_inf(equil, prod);
	// Monte Carlo Integrate:  Integral( x1^2*..*x10^2 * exp( -x1^2 - ...- x10^2)) / Integral ( exp( -x1^2 - ...- x10^2) ) both  from -inf to inf 
	monteCarlo_inf_10D(equil, prod);


	return 0; 
}
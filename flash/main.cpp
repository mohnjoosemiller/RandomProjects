#include "FlashSolver.h" 

/* 
This is the main entry point for the program 
*/ 
int main () 
{
	int nc = 2; 
	double pressure = 1.0 ; 
	double temperature = 300.0; 
	vector<double> zc(nc); 
	zc[0] = 0.2; 
	zc[1] = 1.0 - zc[0];

	// set up solver ( read input, etc ) 
	FlashSolver myFlashSolver(nc); 

	// solve flash at (T,P) 
	int ret = myFlashSolver.solveFlash(pressure, temperature, zc);


	return 0; 
}
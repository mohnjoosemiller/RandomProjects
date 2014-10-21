#include "FlashSolver.h"
#include "SRK_EquationOfState.h"
#include <iostream> 
using namespace std;





FlashSolver::FlashSolver(int nc)
{
	setUpData(nc); 
}

void FlashSolver::setUpData(int nc)
{
	pCompData = new ComponentData(nc);

	cout << "component props set ..." << endl;

}
int FlashSolver::solveFlash(double pres, double temp, vector<double> zc)
{
	int ret = 0; 
	pressure = pres; 
	temperature = temp;
	vector<double> den(2);

	ret = setEOS(SRK_EOS); 

	if ( ret >= 0 ) 
	{
		int ip = 1; //  vap = 0, liq = 1
		
		pEos->computeDensity(zc, pCompData, pres, temp,ip, &den); 

	}
	else 
		return ret;

	return 0;
}

int FlashSolver::setEOS(EOS_T eos )
{
	if ( eos == SRK_EOS )  
		pEos = new SRK_EquationOfState(); 
	else
		return -1; 
	return 0;
}


FlashSolver::~FlashSolver(void)
{
}
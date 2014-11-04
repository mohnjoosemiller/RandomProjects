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

	this->nc = nc;

	cout << "component props set ..." << endl;

}
int FlashSolver::solveFlash(double pres, double temp, vector<double> zc)
{
	int ret = 0; 
	pressure = pres; 
	temperature = temp;
	
	vector<double> k(nc), xi(2*nc), fug(nc);

	double beta = 0.5;

	ret = setEOS(SRK_EOS); 


	pEos->computeFugacity(zc, pCompData, pressure, temperature, 0 , &fug);

	if ( ret >= 0 ) 
	{
		// guess of k-value 
		for ( int i = 0; i < nc; i++ ) 
		{
			double w = pCompData->compOmega[i];
			double tr = temp/(pCompData->compTc[i]);
			double pr = pres/(pCompData->compPc[i]);
			k[i] = (1/pr)*exp(5.37*(1-w)*(1-(1/tr)));
		}

		ret = GetPbubblePdew(pres,temp, zc, k);

		if ( ret >=0) // P is between pbub and pdew
		{	
			cout << "In VLE region, solve flash" << endl;

			ret = RachFordRice(zc, k, &beta, &xi);

			cout << beta << endl;
		}
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

int FlashSolver::RachFordRice(vector<double> zc, vector<double> k, double *beta, vector<double> *xi)
{
	double beta_old = *beta;

	double res = 1, eps = 1e-10;
	int iter = 0; 

	while ( res > eps && iter++ < 100)
	{	
	
		// calculate obj function and derivative
		double objF = 0, objDF = 0; 
		for ( int i = 0; i < nc; i++ ) 
		{
			objF += zc[i]*(k[i]-1)/(1+(beta_old)*(k[i]-1));
			objDF += zc[i]*(k[i]-1)*(k[i]-1)/(1+(beta_old)*(k[i]-1)*(1+(beta_old)*(k[i]-1)));
		}

		(*beta) = beta_old + objF/objDF;

		res = abs((*beta)-beta_old)/abs((*beta));

		beta_old = (*beta);
		cout << iter << "\t" << (*beta) << endl;
	}
	
	if ( iter > 100 || beta < 0)
	{
		cout << "Newton failed to converge or negative beta" << endl;
		cout << "iter = " << iter <<endl;
		return -1;
	}

	// calculate xi & yi
	for ( int i = 0; i < nc; i++)
	{ 
		(*xi)[i] = zc[i]/(1+(*beta)*(k[i]-1));
		(*xi)[i+nc] = k[i]*(*xi)[i];
	}
	return 0;
}

int FlashSolver::GetPbubblePdew(double pres, double temp, vector<double> zc, vector<double> k)
{
	double p_bbl = 0.0; 
	double p_dew = 0.0; 

	vector<double> p_sat(nc); 

	
	

	for (int i = 0; i < nc; i++)
	{
		p_sat[i] = k[i]*pres;

		p_bbl+= zc[i]*p_sat[i];
		p_dew += zc[i]/p_sat[i];
	}
	p_dew = 1./p_dew;

	//cout << p_dew <<"\t"<< p_bbl << endl;

	if ( pres > p_dew && pres < p_bbl )
		return 0;
	return -1; 
}
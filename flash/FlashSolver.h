#include "ComponentData.h"
#include "EquationOfState.h"

enum EOS_T{SRK_EOS};

#pragma once
class FlashSolver
{
public:
	FlashSolver(int nc);
	~FlashSolver(void);
	int setEOS(EOS_T eos); 
	int solveFlash(double p, double t, vector<double> zc);
	int RachFordRice(vector<double> zc, vector<double> k, double *beta, vector<double> *xi);
	int GetPbubblePdew(double pres, double temp, vector<double> zc, vector<double> k);
	void setUpData(int nc);
	double temperature; 
	double pressure;
	ComponentData *pCompData;
	EquationOfState *pEos;
	int nc;
	
};


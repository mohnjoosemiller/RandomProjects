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
	void setUpData(int nc);
	double temperature; 
	double pressure;
	ComponentData *pCompData;
	EquationOfState *pEos;
	
};


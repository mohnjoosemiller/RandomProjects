#include "equationofstate.h"




#pragma once

class SRK_EquationOfState :
	public EquationOfState
{
public:
	SRK_EquationOfState(void);
	~SRK_EquationOfState(void);
	void computeMixParameters(vector<double> zc, ComponentData *p, double Pres, double Temp);
	void computeDensity(vector<double> zc, ComponentData *p,double Pres, double Temp, int phase_id, vector<double> *den);
	void computeFugacity(vector<double> zc, ComponentData *p,double Pres, double Temp, int phase_id, vector<double> *fug);
};


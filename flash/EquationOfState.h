/*
 This is a general equation of state base class, 
 all of the equations of state that are available will implement 
 the required methods slightly differently. 
*/ 
#include "ComponentData.h"
#include <vector> 
#include <iostream> 
using namespace std; 

#define R 83.314

#pragma once
 class EquationOfState
{
public:
	EquationOfState(void);
	~EquationOfState(void);

	/*
	All derived classes should implement these methods
	*/ 
	virtual void computeMixParameters(vector<double> zc, ComponentData *p,double Pres, double Temp)=0;
	virtual void computeDensity(vector<double> zc, ComponentData *p, double Pres, double Temp, int phase_id, vector<double> *den)=0;
	virtual void computeFugacity(vector<double> zc, ComponentData *p, double Pres, double Temp, int phase_id, vector<double> *fug)=0;
	void computeAverageMW(vector<double> zc, ComponentData *p)
	{
			AMW = 0.0; 
			for ( int i = 0; i < p->NC; i++) AMW+=zc[i]*(p->compMW[i]); 
	}; 
	double EOS_B_mix;
	double EOS_A_mix;
	double c0, c1, c2, c3;
	double AMW;

	void solveCubicEOS(double c0, double c1, double c2, double c3, int findMinRoot, double* Zroot );

	ComponentData *pCompData;
};


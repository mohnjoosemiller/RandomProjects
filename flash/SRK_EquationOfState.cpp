#include "SRK_EquationOfState.h"


SRK_EquationOfState::SRK_EquationOfState(void)
{
	cout << "srk set up" << endl;
}


void SRK_EquationOfState::computeMixParameters(vector<double> zc, ComponentData *pComp, double Pres, double Temp)
{
	int nc = pComp->NC; 

// compute average mw 
	computeAverageMW(zc, pComp); 

// compute A and B for mixture, using srk equation
	// first get a, alpha=f(T) and b for pure components
	vector <double> a(nc), alpha(nc), b(nc); 
	double am=0, bm=0, Am=0, Bm=0;
	for ( int i = 0; i < nc; i++ ) 
	{
		b[i] = 0.08664*R*(pComp->compTc[i]) / (pComp->compPc[i]);
		//simple linear mixing rule for bm
		bm += zc[i]*b[i]; 

		a[i] = 0.427*R*R*(pComp->compTc[i])*(pComp->compTc[i]) / (pComp->compPc[i]);
		alpha[i] = (1 + ( 0.48508 + 1.55171*(pComp->compOmega[i])-
			            0.15613*(pComp->compOmega[i])*(pComp->compOmega[i]))*
						(1-sqrt(Temp/(pComp->compTc[i]))));
		alpha[i] *= alpha[i];

		pComp->compai[i] = alpha[i]*a[i];
	}

	for ( int i = 0; i < nc; i++)
		for ( int j = 0; j < nc; j++)
		{
			am += zc[i]*zc[j]*sqrt( alpha[i]*a[i]*alpha[j]*a[j] );
		}

	// Reduced Am and Bm
	Am = am*Pres / ( R*R*Temp*Temp);
	Bm = bm*Pres / (R*Temp); 
	EOS_A_mix = Am; 
	EOS_B_mix = Bm;
	

	// Compute EOS coefficients for cubic eos 
	c0 =  1.0;
	c1 = -1.0;
	c2 = (Am-Bm-Bm*Bm); 
	c3 = -Am*Bm; 
}

void SRK_EquationOfState::computeDensity(vector<double> zc, ComponentData *p, double Pres, double Temp, int phase_id, vector<double> *den)
{
	double Zreturn = 0.0; 

	// update mix parameters 
	computeMixParameters(zc, p, Pres, Temp);

	// solve cubic
	solveCubicEOS(c0,c1,c2,c3, phase_id, &Zreturn);

	double density  = Pres/(Zreturn*R*Temp);

	(*den)[phase_id] = density*AMW;


};


void SRK_EquationOfState::computeFugacity(vector<double> zc, ComponentData *p, double Pres, double Temp, int phase_id, vector<double> *fug)
{

	double Zreturn = 0.0; 

	// update mix parameters 
	computeMixParameters(zc, p, Pres, Temp);
	
	// solve cubic
	solveCubicEOS(c0,c1,c2,c3, phase_id, &Zreturn);

	// get dimensioned mix parameters 
	double bm = EOS_B_mix*R*Temp/Pres;
	double am = EOS_A_mix*R*R*Temp*Temp/Pres;

	// temporary storage
	double term1=0; double term2 = 0; 
	vector<double> lnphi(p->NC);
	
	// compute mixture fugacity coefficients
	for ( int i = 0; i <p->NC; i++)
	{
		// compute "a sum" term 
		double sum_a = 0.0; 
		for ( int j = 0; j < p->NC; j++ )
			sum_a +=zc[j]*sqrt(p->compai[i] * p->compai[j]); 

		term1 = ((p->compbi[i])/bm)*(Zreturn-1)-log(Zreturn-EOS_B_mix);
		term2 = -(EOS_A_mix/EOS_B_mix )*(2*sum_a/am -(p->compbi[i])/bm)*log((Zreturn+EOS_B_mix)/Zreturn); 

		lnphi[i] = term1 + term2;

		// compute fugacity 
		(*fug)[i] = exp(lnphi[i])*zc[i]*Pres;
	} 

};


SRK_EquationOfState::~SRK_EquationOfState(void)
{
}

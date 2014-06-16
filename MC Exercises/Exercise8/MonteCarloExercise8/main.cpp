/* 
Heath Henley
University of Rhode Island

Simple 1D Lennard Jones system 

*/
#include <stdlib.h>
#include <time.h> 
#include <math.h>
#include <iostream>
using namespace std;
#include "Molecule.h"

// kb in kcal/K*molecule
#define K 0.0019872041 

double getEnergy(Molecule &mol)
{
	double r6 = pow((1/mol.position_r),6);
	double r12 = pow( r6, 2);

	return 4* (  r12 - r6  );
}

void moveMolecule(Molecule &mol, double delta, int step, int &accepted)
{
	double random = rand()/RAND_MAX;

	double x_r = mol.position_r;
	double x_new_r = mol.position_r + (random-0.5)*delta;
	//calculate energy - old
	double U = getEnergy(mol);
	//calculate energy - new 
	mol.position_r = x_new_r;
	double U_new = getEnergy(mol);

	if ( x_new_r < 3.0)
	{
		if ( U_new < U ) 
			accepted++;
		else if ( (double)rand()/RAND_MAX < exp ( -mol.beta_r*(U_new-U)) )
			accepted++;
		else // reject 
			mol.position_r = x_r;
	}
	else // reject 
		mol.position_r = x_r;

}


void main () 
{
	srand(time(NULL));
	int Ntrials = 4, prod = 5e6 , equil = 1e6;
	int Nsteps = prod + equil;
	double T = 473.15;
	double beta = 1.0/(K*T), eps = 10, sig = 1;
	double delta = 1;
	int accepted = 0; 
	double Usum = 0.0;

	Molecule mol( eps, sig, sig, beta);
	
	double Usum_trial = 0.0;
	for ( int i = 0; i < Ntrials; i++)
	{
		Usum = 0.0; accepted = 0; delta = 1;
		mol.initialize(eps, sig, sig, beta);

		cout << endl << "Set " << i+1 << " of " << Ntrials << endl; 

		for ( int j = 0; j < Nsteps; j++)
		{
			// try to move molecule
			moveMolecule(mol, delta, j, accepted);

			// get energy
			Usum  += getEnergy(mol);

			if ( j < equil ) 
			{
				if ( j % 100 == 0 ) 
				{
					if ( ((double)accepted / j) < 0.5 )
						delta = delta * 0.95 ;
					else 
						delta = delta * 1.05 ;
				}
			}
			if ( j == equil) 
			{
				cout << "Switch to prod ... " << endl;
				Usum = 0.0;
				accepted = 0;
				cout << "Delta = " << delta << endl;

			}

		}
		Usum_trial += (Usum/prod);
		cout << "Average Dimensionless Energy for set " << i <<"\t"<< Usum/prod << endl;
		cout << "Accepted " << accepted << " moves out of " << prod << ", "<< (double)accepted/prod *100 << " %"<< endl;
	}
	cout << endl<<"Overall Average Energy (dimensionless) = " << Usum_trial/Ntrials << endl;
}
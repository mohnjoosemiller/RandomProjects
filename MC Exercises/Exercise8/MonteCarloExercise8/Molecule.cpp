#include "Molecule.h"


Molecule::Molecule( double e, double s, double p, double beta)
{
	initialize(e,s,p,beta);	
}
void Molecule::initialize( double e , double s, double p, double beta)
{
	eps = e; 
	sigma = s;
	position = p;
	position_r = position/sigma;
	beta_r = beta*eps;
}
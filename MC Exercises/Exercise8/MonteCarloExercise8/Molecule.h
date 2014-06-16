#pragma once
class Molecule
{
public:
	public:
	Molecule(double eps, double sigma, double position, double beta) ;
	void initialize(double eps, double sigma, double position, double beta) ;
	double eps; 
	double sigma; 
	double position;
	double position_r;
	double beta_r;
};


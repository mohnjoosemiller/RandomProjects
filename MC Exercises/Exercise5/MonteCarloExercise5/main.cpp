#include<stdlib.h>
#include<time.h>
#include<iostream>
#include<fstream>
#include<vector> 
using namespace std; 
/* 
Heath Henley 
University of Rhode Island
2014

Two dimensional non interacting lattice gas simulation

*/ 
void getUserInput(int &dx, int &dy, int &N, int &M, int &Nsteps) 
{
// Get all input from user 
	// number of lattice sites in each direction 
	cout << "Enter the number of lattice sites in the x direction: " << endl;
	cin >> dx;
	cout << "Enter the number of lattice site in the y direction: " << endl;
	cin >> dy;
	// total number of lattice sites
	N = dx*dy;
	cout << "\nTotal lattice sites = " << N << endl;
	// number of absorbed molecules
	cout << "Enter the number of absorbed molecules: " << endl;
	cin >> M ;
	cout << "Enter the number of steps: " << endl;
	cin >> Nsteps ;
}

void initLattice(int M, int dx, int dy, vector< vector<int> > &lattice, 
										vector<int> &mols_x,
										vector<int> &mols_y, 
										vector<int> &dx_mol, 
										vector<int> &dy_mol) 
{
	for ( int i = 0; i < dx; i++)
		for (int j = 0; j < dy; j++)
			lattice[i][j] = 0; 

	for ( int i = 0; i < M; i++)
	{
		mols_x[i] = mols_y[i] = 0;
		dx_mol[i] = dy_mol[i] = 0;
	}

	
	int mol_placed = 0; 
	int x_try, y_try;
	while(mol_placed < M)
	{
		x_try = rand() % dx;
		y_try = rand() % dy; 
		// random lattice position is open, so place molecule there
		if( lattice[x_try][y_try] == 0) 
		{
				lattice[x_try][y_try] = 1;
				mols_x[mol_placed] = x_try; 
				mols_y[mol_placed] = y_try;
				mol_placed++;
		}

	}
}

/* print lattice */ 
void printLattice(int dx, int dy, vector< vector<int> > lattice)
{
	for( int i = 0; i < dx; i++) 
	{ 
		for ( int j = 0; j < dy; j++)
			cout << lattice[i][j] << "\t"; 
		cout << endl;
	}
	cout << endl;
}

/* attempt mc move */ 
void performMove(int M, int dx, int dy, vector< vector<int> > &lattice, 
										vector<int> &mols_x,
										vector<int> &mols_y, 
										vector<int> &dx_mol, 
										vector<int> &dy_mol, 
										int & count_accepted) 
{
	// pick random molecule 
	int choose_mol = rand() % M; 
	// pick random direction { -x, +x, -y, +y } 
	int choose_direction = rand() % 4;

	// try move 
	int old_pos_x, new_pos_x;
	int old_pos_y, new_pos_y;

	// save old positions - initalize new position variables
	old_pos_x = mols_x[choose_mol];  
	new_pos_x = old_pos_x; 
	old_pos_y = mols_y[choose_mol];  
	new_pos_y = old_pos_y; 

	if ( choose_direction == 0 ) // -x direction
	{
		new_pos_x = old_pos_x-1; 
		if ( new_pos_x < 0 ) new_pos_x = dx-1;
	}
	else if ( choose_direction == 1 ) // +x direction
	{
		new_pos_x = old_pos_x+1; 
		if ( new_pos_x > dx-1 ) new_pos_x = 0;
	}
	else if ( choose_direction == 2 ) // -y direction
	{	
		new_pos_y = old_pos_y-1; 
		if ( new_pos_y < 0 ) new_pos_y = dy-1;
	}
	else if ( choose_direction == 3 ) // +y direction
	{
		new_pos_y = old_pos_y+1; 
		if ( new_pos_y > dy-1 ) new_pos_y = 0;
	}

	//cout << "Trying to move molecule at {x, y}: {" << old_pos_x <<", "<<old_pos_y<<"}"<<endl;
	//cout << "to coordinates: {" << new_pos_x <<", "<<new_pos_y<<"}"<<endl;

	// trial position has been found, check to see if it is empty 
	// before accepting move 
	if ( lattice[new_pos_x][new_pos_y] == 0 ) 
	{
		// new lattice site is empty, accept move 
		lattice[old_pos_x][old_pos_y] = 0; 
		lattice[new_pos_x][new_pos_y] = 1; 
		mols_x [ choose_mol ] = new_pos_x; 
		mols_y [ choose_mol ] = new_pos_y; 

		// move is accepted, increment the correct displacement counter
		if ( choose_direction == 0 )
		{
			dx_mol[choose_mol] -= 1; 
		}
		else if ( choose_direction ==1)
		{
			dx_mol[choose_mol] += 1; 
		}
		else if ( choose_direction ==2)
		{
			dy_mol[choose_mol] -= 1; 
		}
		else if ( choose_direction ==3)
		{
			dy_mol[choose_mol] += 1; 
		}
		count_accepted ++; 
	}

}

void getMeanDisplacement(int M,  vector<int> dx_mol, vector<int> dy_mol, double &r2)
{
	r2 = 0.0;
	for ( int i = 0 ; i < M; i++) 
	{
		r2 += dx_mol[i]*dx_mol[i] + dy_mol[i]*dy_mol[i];
	}
	r2 = (r2/M);
}

/* 
Main entry point for program
*/ 
int main()
{
// constants and required variables
	int dx, dy, N, M, Nsteps, accepted = 0, Ntrials;
	double r2 = 0.0; 
	bool print_the_lattice = false;

// ouput file 
	fstream fout("r_vs_nstep.txt",ios::out); 

// initialize random seed
	srand(time(NULL));

// Get all input from user 
	//getUserInput( dx, dy, N, M, Nsteps);
	dx = 10; dy  = 10; N = dx*dy; M = 90, Nsteps = 100, Ntrials = 50;

// Coverage = M/N; 
	double theta = (double)M/N;
	cout << "Initial coverage = " << theta << endl;
	
// 2D vector representing 2D lattice
	vector< vector<int> > lattice(dx,dy);

// 1D vectors containing each molecules current position coords
	vector<int> mols_x(M), mols_y(M);  

// 1D vectors constaing overall displacements in each direction
	vector<int> dx_mol(M), dy_mol(M);  

while (Nsteps < 10000)
{
	double ens_sum_r2 = 0.0;
	double ens_sum_r2_sqr = 0.0;

	cout << "Nsteps = " << Nsteps << endl;

	for ( int j = 0; j < Ntrials; j++)
	{
			cout << "\tTrial " <<j<< " of " << Ntrials <<endl;

		// place M molecules randomly on lattice
			initLattice(M, dx, dy, lattice, mols_x, mols_y, dx_mol, dy_mol);
	
		// print initial lattice to debug 
			if (print_the_lattice) printLattice(dx, dy, lattice); 

		// main mc loop 
			double sum_r2 = 0.0;
			for ( int i = 0; i < Nsteps; i++)
			{
				// attempt move 
				performMove(M, dx, dy, lattice, mols_x, mols_y, dx_mol, dy_mol, accepted);
		
				// print to debug
				if (print_the_lattice) printLattice(dx, dy, lattice);

				// get  mean displacement / molecules
				getMeanDisplacement(M,  dx_mol, dy_mol, r2); 
				sum_r2 += r2; 

			}

			double r2_av = sum_r2/Nsteps;
			ens_sum_r2 += r2_av;
			ens_sum_r2_sqr += r2_av*r2_av;
	} // for j ... Ntrials

    double ens_av_r2 = ens_sum_r2/Ntrials;
	double ens_av_r2_sqr = ens_sum_r2_sqr/Ntrials;
	double ens_std = sqrt(ens_av_r2_sqr-ens_av_r2)/(Ntrials-1);
	//cout << "Mean square displacement = " << ens_av_r2 << ",\tStan. Dev = "<<ens_std<< ",\tSteps = "<<Nsteps <<endl;

	fout <<  ens_av_r2 << "\t "<<ens_std<< "\t"<<Nsteps <<endl;

	Nsteps += 200;
}//nstep loop
	return 0 ; 
}
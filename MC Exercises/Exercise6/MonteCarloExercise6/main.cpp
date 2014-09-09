/* 
Heath Henley 
University of Rhode Island
2014

Two dimensional near neighbor interacting lattice gas simulation

*/ 
#include<stdlib.h>
#include<string>
#include<time.h>
#include<iostream>
#include<fstream>
#include<vector> 
using namespace std; 

// k 
const double k = 1.3806488e-23;

// interaction energy
const double EPS =  -(1.3806488e-23*1000);

// global, temp and beta
double T = 273.0;
double BETA = 1/(k*T);
/*
Add frame to movie 
*/ 
void makeMovie(int M, vector<int> mols_x, vector<int> mols_y, fstream &movie)
{
	movie << M << endl<<endl;
	for ( int i = 0; i < M ; i++)
		movie << "H\t"<<mols_x[i] << "\t" << mols_y[i] << "\t" << 0.0<< endl;
}


/* 
Calculate energy of a given configuration
*/
void getEnergy(int dx, int dy, int M, double &E, vector<int> mol_x, vector<int> mol_y)
{
	E = 0.0;
	int x_ref,x;
	int y_ref,y;

	for(int i = 0; i < M-1; i++)
	{
		// position of reference molecule
		x_ref = mol_x[i];
		y_ref = mol_y[i];
		
		// loop through other molecules in system
		for(int j = i+1; j < M; j++)
		{
			x = mol_x[j];
			y = mol_y[j];

			int rx = abs(x-x_ref);
			int ry = abs(y-y_ref);

			// both away from boundary
			if (rx < (dx-1) && ry < (dy-1) ) 
			{
				if ( rx == 1 && ry == 0 || ry == 1 && rx == 0 ) E += EPS;
			} 
			else if ( rx-(dx-1) == 0 && y_ref == y )// both on x boundary
			{
				E += EPS;
			} 
			else if (ry-(dy-1) == 0 && x_ref == x )// both on y boundary
			{
				E += EPS;
			}
		}
	}
}		

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
										int & count_accepted,
										double &E) 
{
	// local energy variables
	double Eold, Enew;
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

	// calculate energy of current configuration
	getEnergy(dx, dy, M, Eold, mols_x, mols_y);

	E = Eold;

	// trial position has been found, check to see if it is empty 
	// before accepting move 
	if ( lattice[new_pos_x][new_pos_y] == 0 ) 
	{

		lattice[old_pos_x][old_pos_y] = 0; 
		lattice[new_pos_x][new_pos_y] = 1; 
		mols_x [ choose_mol ] = new_pos_x; 
		mols_y [ choose_mol ] = new_pos_y; 

		// calculate energy of new configuration
		getEnergy(dx, dy, M, Enew, mols_x, mols_y);

		double test_rand = (double)rand()/RAND_MAX;

		// accept move if energy is lowered or exp(-beta*dE) > rand()
		if ( (Enew < Eold) || exp(-BETA*(Enew-Eold)) > test_rand )
		{
			E = Enew;
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
		else // reject move, reset configuration
		{
			lattice[old_pos_x][old_pos_y] = 1; 
			lattice[new_pos_x][new_pos_y] = 0; 
			mols_x [ choose_mol ] = old_pos_x; 
			mols_y [ choose_mol ] = old_pos_y;
		}
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
	int dx, dy, N, M, Nsteps, accepted = 0, Ntrials, equil, prod, movieFreq;
	double r2 = 0.0, E = 0.0;  
	bool print_the_lattice =false;

// ouput file 
	fstream fout("out.txt",ios::out); 

// movie file 
	fstream movie;
	string movie_file;

// initialize random seed
	srand(time(NULL));

// Get all input from user 
	//getUserInput( dx, dy, N, M, Nsteps);
	equil = 2000; 
	prod  = 10000;
	dx = 10; dy  = 10;
	N = dx*dy; 
	M = 50; 
	Nsteps = equil+prod; 
	Ntrials = 4;
	movieFreq = 10;

// Coverage = M/N; 
	double theta = (double)M/N;
	cout << "Initial coverage = " << theta << endl;
	
// 2D vector representing 2D lattice
	vector< vector<int> > lattice(dx,dy);

// 1D vectors containing each molecules current position coords
	vector<int> mols_x(M), mols_y(M);  

// 1D vectors constaing overall displacements in each direction
	vector<int> dx_mol(M), dy_mol(M);  


//while (T < 400.0)
//{
	double ens_sum_r2 = 0.0;
	double ens_sum_r2_sqr = 0.0;
	double ens_sum_E = 0.0;
	double ens_sum_E2 = 0.0; 

	for ( int j = 0; j < Ntrials; j++)
	{
			movie_file = "movie_set" + to_string((long long int)j) + ".xyz" ;
			movie.open(movie_file, ios::out);

			cout <<endl << "Begin Trial " <<j<< " of " << Ntrials <<endl;
			cout << "Nstep \t Energy/part (K)" <<endl;

		// place M molecules randomly on lattice
			initLattice(M, dx, dy, lattice, mols_x, mols_y, dx_mol, dy_mol);
	
		// print initial lattice to debug 
			if (print_the_lattice) printLattice(dx, dy, lattice); 

		// main mc loop 
			double sum_r2 = 0.0, sum_E = 0.0;
			for ( int i = 0; i < Nsteps; i++)
			{

				if ( i == equil)
				{
					cout << "Switch to production mode ... "  << endl;
					sum_r2 = 0.0; 
					sum_E = 0.0; 
				}

				// attempt move 
				performMove(M, dx, dy, lattice, mols_x, mols_y, dx_mol, dy_mol, accepted,E);

				// print to debug
				if (print_the_lattice) printLattice(dx, dy, lattice);

				// get  mean displacement / molecules
				getMeanDisplacement(M,  dx_mol, dy_mol, r2); 
				sum_r2 += r2;

				// get and sum energy
				getEnergy( dx, dy, M, E, mols_x, mols_y);
				sum_E += E/k;

				if ( i % 2000 == 0 ) 
				{
					cout << i << "\t" << E/(M*k) << endl;
				}
				if ( i % movieFreq == 0 ) 
					makeMovie(M, mols_x, mols_y, movie);

			}
			movie.close();
			// mean square displacement
			double r2_av = sum_r2/prod;
			ens_sum_r2 += r2_av;
			ens_sum_r2_sqr += r2_av*r2_av;
			//energy
			double E_av = sum_E/(prod*M);
			ens_sum_E +=  E_av;
			ens_sum_E2 += E_av*E_av;


	} // for j ... Ntrials
	// mean square displacement 
    double ens_av_r2 = (ens_sum_r2/Ntrials);
	double ens_av_r2_sqr = ens_sum_r2_sqr/Ntrials;
	double ens_std = sqrt(ens_av_r2_sqr-ens_av_r2*ens_av_r2)/(Ntrials-1);
	// energy 
	double E_tot_av = ens_sum_E/Ntrials;
	double E2_tot_av = ens_sum_E2/Ntrials;
	double std_E = sqrt(E2_tot_av-E_tot_av*E_tot_av)/(Ntrials-1);

	
	cout <<endl<< "Final Results for " << Ntrials << " trials:" << endl;
	cout << Nsteps << "\t"<<theta<< "\t" << T <<"\t"<< ens_av_r2 << "\t "<<ens_std<< "\t " << E_tot_av << "\t " << std_E << endl;

	fout << Nsteps << "\t"<<theta << "\t" <<T <<"\t"<< ens_av_r2 << "\t "<<ens_std<< "\t " << E_tot_av << "\t " << std_E << endl;

	//T += 30;
	//BETA = 1.0/(k*T); 
// }//T loop

	return 0 ; 
}
#include "ComponentData.h"

#define NC_MAX 3

// data is not correct yet, but dummied in to debug
// critical temperature in K 
double tc_lib[NC_MAX] = {
	0.0, 
	647.15 , 
	304.1 };

// critical pressure in bar 
double pc_lib[NC_MAX] = {
	0.0, 
	200, 
	73 };

// critical volume in cm3/mol 
double vc_lib[NC_MAX] = {
	0.0, 
	18, 
	34 };

// critical compressiblity  
double zc_lib[NC_MAX] = {
	0.0, 
	0.27, 
	0.27 };

// accentric factor 
double omega_lib[NC_MAX] = {
	0.0, 
	0.244, 
	0.244 };

// molecular weight
double mw_lib[NC_MAX] = {
	0.0, 
	18.0145, 
	44.01 };




ComponentData::ComponentData(int nc)
{
	// read data ( this should be input from a file eventually 
	//int nc = 2; 
	vector<int> ids(nc);
	ids[0] = 1; ids[1] =2; 

	// set nc 
	NC = nc; 

	// initialize arrays
	initializeDataMembers();

	// set component ids
	for ( int i = 0; i < nc ; i++ ) compID[i] = ids[i];

	// set property data
	setPropertyData();
}

void ComponentData::initializeDataMembers()
{
	compID.resize(NC);
	compTc.resize(NC); 
	compPc.resize(NC);
	compVc.resize(NC);
	compZc.resize(NC); 
	compMW.resize(NC);
	compbi.resize(NC); 
	compOmega.resize(NC); 
}

ComponentData::~ComponentData(void)
{
}

void ComponentData::setPropertyData()
{
	for ( int i = 0; i < NC; i++)
	{
		compTc[i] = tc_lib[compID[i]];
		compPc[i] = pc_lib[compID[i]];
		compVc[i] = vc_lib[compID[i]];
		compZc[i] = zc_lib[compID[i]];
		compOmega[i] = omega_lib[compID[i]];
		compMW[i] = mw_lib[compID[i]];
	}
}




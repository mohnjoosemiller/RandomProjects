#include<vector>
using namespace std;

#pragma once
class ComponentData
{
public:
	ComponentData(int nc);
	~ComponentData(void);
	vector<int> compID;
	vector<double> compTc;
	vector<double> compPc;
	vector<double> compVc;
	vector<double> compZc;
	vector<double> compOmega;
	vector<double> compMW;
	vector<double> compbi;
	vector<double> compai;

	int NC;
private:
	void initializeDataMembers();
	void setPropertyData();
};


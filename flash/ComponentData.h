#include<vector>
using namespace std;

#pragma once
class ComponentData
{
public:
	ComponentData(int nc);
	~ComponentData(void);
	vector<int> compID;
	vector<int> compTc;
	vector<int> compPc;
	vector<int> compVc;
	vector<int> compZc;
	vector<int> compOmega;
	vector<int> compMW;
	vector<int> compbi;

	int NC;
private:
	void initializeDataMembers();
	void setPropertyData();
};


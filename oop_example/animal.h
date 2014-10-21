#pragma once
#include <iostream> 
using namespace std; 

class animal
{
public:
	animal(void);
	~animal(void);
	virtual void run() = 0;
	virtual void speak()= 0;
	virtual void move() = 0;
};


#pragma once
#include "animal.h"
class Dog :
	public animal
{
public:
	Dog(void);
	~Dog(void);
	void run();
	void move(); 
	void speak(); 
};


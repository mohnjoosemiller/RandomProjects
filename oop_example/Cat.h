#pragma once
#include "animal.h"
class Cat :
	public animal
{
public:
	Cat(void);
	~Cat(void);
	void run();
	void move(); 
	void speak(); 
};


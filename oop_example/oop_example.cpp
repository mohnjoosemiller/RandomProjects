// oop_example.cpp : Defines the entry point for the console application.
//
#include<stdlib.h> 
#include<iostream>
using namespace std; 

#include "animal.h"
#include "Dog.h"
#include "Cat.h"

enum ANIMALS { CAT, DOG };


int main(int argc, char* argv[])
{
	animal *myAnimal; 

    ANIMALS ANIMAL_T = CAT;

	if ( ANIMAL_T == DOG) 
		myAnimal = new Dog(); 
	else if (ANIMAL_T == CAT)
		myAnimal = new Cat();
	else
		myAnimal = NULL;


	if ( myAnimal != NULL)
	{ 
		myAnimal->move();
		myAnimal->speak();
	}
	return 0;
}


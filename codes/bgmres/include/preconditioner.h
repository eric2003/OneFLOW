#ifndef PRECONDITIONERCLASS
#define PRECONDITIONERCLASS
#include "systemSolver.h"
#include "UCom.h"
class Solution;
class UCom;

class Preconditioner
{

public:
	Preconditioner(int number= Rank.NUMBER);              //< Default constructor for the class
	Preconditioner(const Preconditioner& oldCopy);  //< Constructor for making a copy/duplicate
	~Preconditioner();                              //< Destructor for the class

	Solution solve(const Solution &vector);    //< Method to solve the
	Solution solve2(const Solution &vector);																					 //< system associated with
																						 //< the preconditioner.

	
	/**
		 Method to set the number of elements to use for the length of the approximation.

		 @param number The number of elements.
		 @return N/A
	 */
	void setN(int number) 
	{
		N = number;
	}

	/**
		 Method to get the number of elements that are used for the approximation.

		 @return The number of grid points used in the approximation.
	 */
	int getN() const
	{
		return(N);
	}

	/**
		 Method to get the value of the preconditioner's vector for a
		 given row.

		 @param row Row in the vector you want to access.
		 @return The value of the vector for the given row.
	 */
	double getValue(int row,int col) const
	{
		return(vector[row][col]);
	}


protected:


private:

	int N;           //< The number of grid points associated with the approximation.
	double **vector; //< The vector that has the reciprocol of the diagonal entries of the operator.
};




#endif

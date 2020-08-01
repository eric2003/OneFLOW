#include "systemSolver.h"
#pragma once
class Solution;

class Poisson
{
public:
	friend class UINsInvterm;
public:
	Poisson(int number= Rank.NUMBER, int ranknumber = Rank.RANKNUMBER);        //< Default constructor for the Poisson Class.
	Poisson(const Poisson& oldCopy);   //< Copy constructor for the Poisson Class.
	~Poisson();                        //< Destructor for the Poisson Class.

	// Basic algebraic operators associated with the linearization of the operator.
	double& operator()(int row);     //< The value of the linearization for the operator at a given row and column.
	Solution operator*(class Solution vector);  //< The linearized operator acting on a given Solution.

	/**
		 Method to get the number of elements that are in the approximation.

		 @return The number of elements in the grid.
	 */
	int  getN() const
	{
		return(N);
	}
	 /**
	     Method to get the rank of the matrix that are in the matrices
	 */
	int getL() const
	{
		return(L);
	}

	/**
		 Method to get one of the elements from the first derivative matrix.

		 @param row The row number in the matrix
		 @param col The column number in the matrix.
		 @return The value within the matrix at the given row and column.
	 */
	double geta(int row) const
	{
		return(A[row]);
	}

protected:
	/*void produceA(double* deriv, int* derivrow, int* derivcol);*/
	void coeMatrix(double* deriv, int* derivrow, int* derivcol);

public:
	int N;        //< The number of grid points in the approximation.
	int L;        //< The rank of the matrix
	double *A;    //< The linear operator matrix
	int *IA;      //< The row number of the matrix
	int *JA;      //< The column number of the matrix
};


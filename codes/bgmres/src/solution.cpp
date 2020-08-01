#include <iostream>
#include <cmath>
#include <vector>
#include <stdio.h>
#include "solution.h"
#include "poisson.h"
#include "util.h"

using namespace std;
/** ************************************************************************
 * Base constructor  for the Solution class. 
 * @param size The length of the vector used in the approximation (optional)
 * ************************************************************************ */
Solution::Solution(int size)
{
	// Set the size of the vector, allocate the space, and zero out the
	// approximation.
	setN(size);
	solution = ArrayUtils<double>::twotensor(size,Rank.COLNUMBER);  // allocate the space. 

}

/** ************************************************************************
 * Copy constructor  for the Solution class. 
 * @overload
 * @param oldCopy The Solution class member to make a copy of.
 * ************************************************************************ */
Solution::Solution(const Solution& oldCopy)
{
	// Make a copy of the Solution that is passed to me.
	// Set the size of the vector, allocate the space, and then
	// copy the values over.
	int size = oldCopy.getN();
	setN(size);
	solution = ArrayUtils<double>::twotensor(size,Rank.COLNUMBER);
	for(size = size-1;size>=0;size--)
	{
		for(int col=Rank.COLNUMBER-1;col>=0;col--)
		{
			setEntry(oldCopy.getEntry(size,col),size,col);
		}
	}
}

/** ************************************************************************
 *	Destructor for the Solution class. 
 *  ************************************************************************ */
Solution::~Solution()
{
	// delete the approximation.
	if(solution)
	ArrayUtils<double>::deltwotensor(solution);
	solution = NULL;
}

/** ************************************************************************
 * The parenthesis operator for the Solution class.
 * 
 * Returns the value of the approximation for the indicated row.
 *
 * @param row The row number to use.
 * @return a double precision value, solution[row]
 * ************************************************************************ */
double& Solution::operator()(int row,int col)
{
	return(solution[row][col]);
}

/** ************************************************************************
 * The equals operator for the Solution class.
 * 
 * Returns a new Solution object that is identical to the Solution
 * object passed to it.
 *
 * @param vector The Solution argument to copy
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution Solution::operator=(const Solution& vector)
{
	int row;
	int col;
	int N = getN();
	if(this != &vector)
		{
			for(row=N-1;row>=0;row--)
				for(col=Rank.COLNUMBER-1;col>=0;col--)
					{
						this->setEntry(vector.getEntry(row,col),row,col);
					}
		}
	return(*this);
}

/** ************************************************************************
 * The equals operator for the Solution class.
 * 
 * Returns a new Solution object in which every entry in the vector is
 * the same as the double precision number passed to it.
 *
 * @overload
 * @param value The value to copy into the vector.
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution Solution::operator=(const double& value)
{
	int row;
	int col;
	int N = getN();
	for(row=N-1;row>=0;row--)
		for(col=Rank.COLNUMBER-1;col>=0;col--)
			{
				this->setEntry(value,row,col);
			}

	return(*this);
}

/** ************************************************************************
 * The summation operator for the Solution class.
 * 
 * Returns a new Solution object which is the sum of two Solution
 * objects.
 *
 * @param vector The Solution object to add to the current Solution object.
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution Solution::operator+(const Solution& vector)
{
	int N = vector.getN();
	Solution result(N);
	int row;
	int col;
	
	for(row=N-1;row>=0;row--)
		for(col=N-1;col>=0;col--)
			{
				result.setEntry(this->getEntry(row,col)+vector.getEntry(row,col),row,col);
			}

	return(result);

}

/** ************************************************************************
 * The subtraction operator for the Solution class.
 * 
 * Returns a new Solution object which is the difference of two Solution
 * objects.
 *
 * @param vector The Solution object to subtract from the current Solution object.
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution Solution::operator-(const Solution& vector)
{
	int N = vector.getN();
	Solution result(N);
	int row;
	int col;
	for(row=N-1;row>=0;row--)
		for(col=Rank.COLNUMBER-1;col>=0;col--)
			{
				result.setEntry(this->getEntry(row,col)-vector.getEntry(row,col),row,col);
			}

	return(result);

}

/** ************************************************************************
 * The scalar multiplication operator for the Solution class.
 * 
 * Returns a new Solution object which is the scalar product of an
 * object from the Solution class and a single double precision
 * number.
 *
 * @param value Scalar value to multiply every entry in the current vector.
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution Solution::operator*(const double& value)
{
	int N = getN();
	Solution result(N);
	int row;
	int col;
	for(row=N-1;row>=0;row--)
		for(col=Rank.COLNUMBER-1;col>=0;col--)
			{
				result.setEntry(this->getEntry(row,col)*value,row,col);
			}

	return(result);
}

/** ************************************************************************
 * The dot product operation for the Solution class.
 * 
 * Returns the dot product of this object and another object in the
 * Solution class.
 *
 * @overload
 * @param vector The other solution object to use for the dot product.
 * @return the dot product.
 * ************************************************************************ */
double Solution::operator*(const Solution& vector)
{
	int N = vector.getN();
	double dot = 0.0;
	int row;
	int col;
	for(row=N-1;row>=0;row--)
		for(col=Rank.COLNUMBER-1;col>=0;col--)
			{
				dot += getEntry(row,col)*vector.getEntry(row,col);
			}
	return(dot);
}

/** ************************************************************************
 * The scalar multiplication operator for "*=" for the Solution class.
 * 
 * Returns a new Solution object which is the scalar product of an
 * object from the Solution class and a single double precision
 * number.
 *
 * @param value Scalar value to multiply every entry in the current vector.
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution Solution::operator*=(const double& value)
{
	int N = getN();
	int row;
	int col;
	for(row=N-1;row>=0;row--)
		for(col=N-1;col>=0;col--)
			{
				this->setEntry(this->getEntry(row,col)*value,row,col);
			}

	return(*this);
}

/** ************************************************************************
 * The subtraction operator for "-=" for the Solution class.
 * 
 * Returns a new Solution object which is found by subtracting each
 * element from the passed Solution object from the current object.
 *
 * @param vector The Solution object to subtract from this object.
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution Solution::operator-=(const Solution& vector)
{
	int N = vector.getN();
	int row;
	int col;
	for(row=N-1;row>=0;row--)
		for(col=N-1;col>=0;col--)
			{
				this->setEntry(this->getEntry(row,col)-vector.getEntry(row,col),row,col);
			}

	return(*this);
}

/** ************************************************************************
 * The addition operator for "+=" for the Solution class.
 * 
 * Returns a new Solution object which is found by adding each
 * element from the passed Solution object to the current object.
 *
 * @param vector The Solution object to add to this object.
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution Solution::operator+=(const Solution& vector)
{
	int N = vector.getN();
	int row;
	int col;
	for(row=N-1;row>=0;row--)
		for(col=N-1;col>=0;col--)
			{
				this->setEntry(this->getEntry(row,col)+vector.getEntry(row,col),row,col);
			}

	return(*this);
}

/** ************************************************************************
 * The scalar multiplication operator for the Solution class.
 * 
 * Returns a new Solution object which is found by multiplying each
 * element from the passed Solution object by the scalar value.
 *
 * @relates Solution
 * @param value The scalar value (Double) to multiply.
 * @param vector The Solution object to multiply the scalar by.
 * @return An object from the Solution class.
 * ************************************************************************ */
Solution operator*(const double &value,class Solution vector)
{
	return(vector*value);
}


/** ************************************************************************
 * The method to find the dot product between two solutions.
 * 
 * Takes two solution vectors and computes their dot product. 
 *
 * @param value The first Solution object to use.
 * @param row The second Solution object to use.
 * @return The scalar dot product.
 * ************************************************************************ */
double Solution::dot(const Solution& v1,const Solution& v2)
{
	int N = v1.getN();
	double dotProduct = 0.0;
	int row;
	int col;
	//std::cout << "dot product" << std::endl;
	for(row=N-1;row>=0;row--)
		for(col=N-1;col>=0;col--)
			{
				dotProduct += v1.getEntry(row,col)*v2.getEntry(row,col);
			}
	return(dotProduct);
}

/** ************************************************************************
 * The method to find the dot product between two solutions.
 * 
 * Takes two solution vectors and computes their dot product. 
 *
 * @overload
 * @param value The first Solution object to use.
 * @param row The second Solution object to use.
 * @return The scalar dot product.
 * ************************************************************************ */
double Solution::dot(Solution* v1,Solution* v2)
{

	double dotProduct = 0.0;
	int N = v1->getN();
	int row;
	int col;
	for(row=N-1;row>=0;row--)
		for(col=N-1;col>=0;col--)
			{
				dotProduct += v1->getEntry(row,col)*v2->getEntry(row,col);
			}
	return(dotProduct);
}

/** ************************************************************************
 * The method to find the norm of a solution.
 * 
 * Takes a Solution object and computes the l2 norm of the solution.
 *
 * @param v1 The  Solution object to use.
 * @return The norm of the Solution object.
 * ************************************************************************ */
double Solution::norm(const Solution& v1)
{
	int N = v1.getN();
	double norm = 0.0;
	int row;
	int col;
	for(row=N-1;row>=0;row--)
		for(col=N-1;col>=0;col--)
			{
				norm += v1.getEntry(row,col)*v1.getEntry(row,col);
			}
	return(sqrt(norm));
}

/** ************************************************************************
 * The method to find the norm of a solution.
 * 
 * Determines the l2 norm of the associated object.
 *
 * @overload
 * @return The norm of the Solution object.
 * ************************************************************************ */
double Solution::norm()
{
	int N = getN();
	double norm = 0.0;
	int row;
	int col;
	for(row=N-1;row>=0;row--)
		for(col=Rank.COLNUMBER-1;col>=0;col--)
			{
				norm += getEntry(row,col)*getEntry(row,col);
			}
	return(sqrt(norm));
}

/** ************************************************************************
* The method to find the norm of a single vector.
* ************************************************************************** */
double Solution::normSingle()
{
	double normSingle = getEntry(0,0)*getEntry(0,0);
	int lupe;
	for(lupe=getN()-1;lupe>0;lupe--)
		{
			normSingle += getEntry(lupe,1)*getEntry(lupe,1);
		} 
	return(sqrt(normSingle));
}

/** ************************************************************************
* The method to find the norm of a single vector.
* ************************************************************************** */
double Solution::normSquare()
{
	double normSquare = 0.0;
	int row;
	int col;
	for(row=Rank.COLNUMBER-1;row>=0;row--)
		for(col=Rank.COLNUMBER-1;col>=0;col--)
			{
				normSquare += getEntry(row,col)*getEntry(row,col);
			}
	return(sqrt(normSquare));
}

/** 
* Matrix Multiplied
*/
double Solution::Multiplied(const Solution& v1,const Solution& v2)
{
	int N = v1.getN();
	double dotProduct = 0.0;
	int row;
	int col;
	//std::cout << "dot product" << std::endl;
	for(row=N-1;row>=0;row--)
		for(col=N-1;col>=0;col--)
			{
				dotProduct += v1.getEntry(row,col)*v2.getEntry(row,col);
			}
	return(dotProduct);
}

/** ************************************************************************
 * The method to perform an axpy procedure.
 * 
 * Subtracts multipler*vector from the current approximation.
 *
 * @param The vector to add
 * @param The scalar multiple to add
 * @return N/A
 * ************************************************************************ */
void Solution::axpy(Solution* vector,
					double multiplier)
{
	int N = vector->getN();
	int row;
	int col;
	for(row=N-1;row>=0;row--)
		for(col=N-1;col>=0;col--)
			{
				setEntry(getEntry(row,col)+multiplier*vector->getEntry(row,col),
						 row,col);
			}
}


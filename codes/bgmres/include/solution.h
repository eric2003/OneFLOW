#ifndef SOLUTIONCLASS
#define SOLUTIONCLASS
#include "systemSolver.h"
#include "poisson.h"
#include "util.h"

class Solution
{

public:
	Solution(int size=Rank.RANKNUMBER);               //< Default constructor for the class
	Solution(const Solution& oldCopy);       //< Constructor for making a copy/duplicate
	~Solution();                             //< Destructor for the class

	// Now define the operators associated with the class.
	double& operator()(int row, int col);                   //< The parenthesis operator for access to data elements
	Solution operator=(const Solution& vector);    //< Assignment operator for copying another Solution
	Solution operator=(const double& value);       //< Assignment operator for assigning a single value to all elements.
	Solution operator+(const Solution& vector);    //< Operator for adding two Solution objects
	Solution operator-(const Solution& vector);    //< Operator for subtracting two Solution objects.
	Solution operator*(const double& value);       //< Operator for scalar multiplication.
	double   operator*(const Solution& vector);    //< Operator for the dot product
	Solution operator*=(const double& value);      //< Operator for scalar multiplication in place.
	Solution operator-=(const Solution& vector);   //< Operator for subtracting another Solution object.
	Solution operator+=(const Solution& vector);   //< Operator for adding another Solution object.


	/** Definition of the dot product of two approximation vectors. */
	static double dot (const Solution& v1,const Solution& v2);
	static double dot (Solution* v1,Solution* v2);

	/** Definition of the l2 norm of an approximation vector. */
	static double norm(const Solution& v1);
	double norm();
	double normSingle();
	double normSquare();
    static double Multiplied (const Solution& v1,const Solution& v2);

	/** Definition of the axpy procedure. */
	void axpy(Solution* vector,double multiplier);
	//static double getcol(const Solution& vector);

	/** ************************************************************************
	 * The method to set the value of the entry in a row of the solution.
	 * 
	 * Sets the value of the indicated row to the value specified.
	 *
	 * @param value The scalar value (double) set the given row to.
	 * @param row The entry in the vector to change.
	 * @return N/A
	 * ************************************************************************ */
	void setEntry(double value,int row, int col)
	{
		solution[row][col] = value;
	}



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
	   Method to get the number of elements used for the approximation.

	   @return The number of elements in the approximation.
	*/
	inline int getN() const
	{
		return(N);
	}

	/**
	   Method to get the value of the approximation at a certain grid point.

	   @param row The grid point where you want the height of the function.
	   @return The approximation at the given grid point.
	*/
	inline double getEntry(int row, int col) const
	{
		return(solution[row][col]);
	}

protected:



private:

	// Define the size of the vector and the vector that will contain
	// the information.
	int N;                      //< The number of grid points.
	double **solution ;    //< The vector that contains the approximation.

};


/**
	 Method to define scalar multiplcation for a solution vector.

	 The method to define how to multiply an approximation on the left using a scalar product.

	 @param value  The scalar to multiply the approximation by
	 @param vector The approximation that is being multiplied by the scalar.
	 @return The result of the scalar multiplication on the approximation.
 */
Solution operator*(const double& value,class Solution vector);

#endif

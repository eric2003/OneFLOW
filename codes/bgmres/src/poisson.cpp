
#include "poisson.h"
#include "solution.h"
#include "util.h"
#include "iostream"
#include "fstream"
#include <cmath>
#include <vector>
#include <stdio.h>
#include <UCom.h>
#include "UINsInvterm.h"
#include "INsInvterm.h"



Poisson::Poisson(int number, int ranknumber)
{
    N  = number;
	L  = ranknumber;
	A  = ArrayUtils<double>::onetensor(number);
	IA = ArrayUtils<int>::onetensor(number);
	JA = ArrayUtils<int>::onetensor(number);
	coeMatrix(A,IA,JA);  //按行存储方法
}


Poisson::Poisson(const Poisson& oldCopy)
{
    N  = oldCopy.getN();
	L  = oldCopy.getL();
	A  = ArrayUtils<double>::onetensor(N);
	IA = ArrayUtils<int>::onetensor(N);
	JA = ArrayUtils<int>::onetensor(N);
}

Poisson::~Poisson()
{
	ArrayUtils<double>::delonetensor(A);
	ArrayUtils<int>::delonetensor(IA);
	ArrayUtils<int>::delonetensor(JA);
}

double& Poisson::operator()(int row)
{
	return(A[row]);
}

/**
 * 对应coeMatrix压缩函数
*/
Solution Poisson::operator*(class Solution vector)
{
	
	int i, j;
	int a, b;
	double data;
	double** tmp = ArrayUtils<double>::twotensor(Rank.RANKNUMBER, Rank.COLNUMBER);
	Solution result;
	// the first and last row just return the same values, so 
	// there is no need to define the results from that row.
	
	for( i = 1; i <= Rank.RANKNUMBER; i++)               
	{
		for(j = 0; j < Rank.COLNUMBER; j++)
	    {
			for(int col = IA[i-1]; col < IA[i]; col++)
			{
				b = JA[col];
				data = A[col];
				tmp[i-1][j] += data * vector(b,j);
				result.setEntry(tmp[i-1][j],i-1,j);
				/*cout << tmp[a-1][j] << endl;*/
			}
		}
	}
	ArrayUtils<double>::deltwotensor(tmp);
	return(result);
}

void Poisson::coeMatrix(double* deriv,int* derivrow,int* derivcol)
{
	for(int m = 0; m <= Rank.RANKNUMBER; ++m)                 
	{
			derivrow[m] = Rank.TempIA[m];
	}
	for (int k = 0; k < Rank.NUMBER; ++k)
	{
		derivcol[k] = Rank.TempJA[k];
		deriv[k] = Rank.TempA[k];
	}
}

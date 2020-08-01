#include "poisson.h"
#include "solution.h"
#include "preconditioner.h"
#include "GMRES.h"
#include "fstream"
#include <iostream>
#include <cmath>
#include <time.h>
#include <stdlib.h>
#include "systemSolver.h"
#include "UCom.h"
#include <UINsInvterm.h>

SolveMRhs bgx;
SolveMRhs::SolveMRhs()
{
	;
}

SolveMRhs::~SolveMRhs()
{
	;
}

SolveMRhs Rank;
void SolveMRhs::Init()
{
	TempA = ArrayUtils<double>::onetensor(Rank.NUMBER);
	TempIA = ArrayUtils<int>::onetensor(Rank.RANKNUMBER+1);
	TempJA = ArrayUtils<int>::onetensor(Rank.NUMBER);
	TempB = ArrayUtils<double>::twotensor(Rank.RANKNUMBER,Rank.COLNUMBER);
	TempX = ArrayUtils<double>::twotensor(Rank.RANKNUMBER,Rank.COLNUMBER);
}
void SolveMRhs::Deallocate()
{
	ArrayUtils<double>::delonetensor(TempA);
	ArrayUtils<int>::delonetensor(TempIA);
	ArrayUtils<int>::delonetensor(TempJA);
	ArrayUtils<double>::deltwotensor(TempB);
	ArrayUtils<double>::deltwotensor(TempX);
}
void SolveMRhs::BGMRES()
{
	clock_t start, finish;
	double time;
	start = clock();
	Poisson* A = new Poisson;   // The operator to invert.
	Solution* x = new Solution(Rank.RANKNUMBER);  // The approximation to calculate.
	Solution* b = new Solution(Rank.RANKNUMBER);  // The forcing function for the r.h.s.
	Solution* residual = new Solution(Rank.RANKNUMBER);
	Preconditioner* pre =
		new Preconditioner(Rank.RANKNUMBER);      // The preconditioner for the system.
	int restart = 0;                    // Number of restarts to allow
	int maxIt = 500;                      // Dimension of the Krylov subspace
	double tol = 1.0E-8;                 // How close to make the approximation.

	/**
	   produce the right-hand sides
	*/
	int i, j;
	for (i = 0; i < Rank.RANKNUMBER; i++)
	{
		for (j = 0; j < Rank.COLNUMBER; j++)
		{
			{
				(*b)(i, j) = Rank.TempB[i][j];
			}
		}
	}
	// Find an approximation to the system!
	int result = GMRES(A, x, b, residual, pre, maxIt, restart, tol);

	// Output the solution
	for (int lupe = 0; lupe < Rank.COLNUMBER; lupe++)
	{
		for (int innerlupe = 0; innerlupe < Rank.RANKNUMBER; innerlupe++)
		{
			Rank.TempX[innerlupe][lupe] = (*x)(innerlupe, lupe);
		}
	}

	//std::cout << "Iterations: " << result << " residual: " << tol << std::endl;
	finish = clock();
	time = (double)(finish - start);    //计算运行时间

	delete A;
	delete x;
	delete b;
	delete residual;
	delete pre;


#define SOLUTION
#ifdef SOLUTION
	/*ofstream file("solution.txt", ios::app);
	int lupe;
	int innerlupe;
	for (lupe = 0; lupe < Rank.RANKNUMBER; ++lupe)
	{
		for (innerlupe = 0; innerlupe < Rank.COLNUMBER; ++innerlupe)
		{
			file << lupe << "," << (*x)(lupe, innerlupe)
				<< endl;
		}
	}
	file << time << endl;
	file.close();*/
	
#endif
}

#include "poisson.h"
#include <memory>
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
	auto A = std::make_unique<Poisson>();   // The operator to invert.
	auto x = std::make_unique<Solution>(Rank.RANKNUMBER);  // The approximation to calculate.
	auto b = std::make_unique<Solution>(Rank.RANKNUMBER);  // The forcing function for the r.h.s.
	auto residual = std::make_unique<Solution>(Rank.RANKNUMBER);
	auto pre = std::make_unique<Preconditioner>(Rank.RANKNUMBER);  // The preconditioner for the system.
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
	int result = GMRES(A.get(), x.get(), b.get(), residual.get(), pre.get(), maxIt, restart, tol);

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
	time = (double)(finish - start);    //Calculate run time

	// unique_ptr members destroy A, x, b, residual, pre.


#define SOLUTION
#ifdef SOLUTION
	/*ofstream file("solution.txt", std::ios::app);
	int lupe;
	int innerlupe;
	for (lupe = 0; lupe < Rank.RANKNUMBER; ++lupe)
	{
		for (innerlupe = 0; innerlupe < Rank.COLNUMBER; ++innerlupe)
		{
			file << lupe << "," << (*x)(lupe, innerlupe)
				<< std::endl;
		}
	}
	file << time << std::endl;
	file.close();*/
	
#endif
}

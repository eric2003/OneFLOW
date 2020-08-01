
#include "preconditioner.h"
#include "solution.h"
#include "util.h"
#include "UCom.h"
#include <cmath>

/** ************************************************************************
 * Base constructor  for the Preconditioner class. 
 *
 *
 * @param size The number of grid points used in the approximation.
 * ************************************************************************ */
Preconditioner::Preconditioner(int number)
{
	setN(number);
}

/** ************************************************************************
 *	Copy constructor  for the Preconditioner class. 
 *
 *	@param oldCopy The Preconditioner class member to make a copy of.
 * ************************************************************************ */
Preconditioner::Preconditioner(const Preconditioner& oldCopy)
{
	setN(oldCopy.getN());
}

/** ************************************************************************
 *	Destructor for the Preconditioner class. 
 *  ************************************************************************ */
Preconditioner::~Preconditioner()
{
	;
	//ArrayUtils<double>::deltwotensor(vector);
}

/** ************************************************************************
 * The method to solve the system of equations associated with the
 * preconditioner.
 * 
 * Returns the value Solution class that is the solution to the
 * preconditioned system.  For the moment we just assume that the
 * preconditioner is the identity matrix. (This needs to be changed!)
 *
 * @param vector The Solution or right hand side of the system.
 * @return A Solution class member that is the solution to the
 *         preconditioned system.
 * ************************************************************************ */
Solution Preconditioner::solve(const Solution &current)
{
	Solution multiplied(current);
	int cId, iFace, lc, rc, offdiag, offdiag2;
	int innerlupe;
	int fId, fn1, n1, n2;
	int fn2 = 0;
	for (cId = 0; cId < getN(); cId++)
	{
		for (innerlupe = 0; innerlupe < Rank.COLNUMBER; innerlupe++)
		{
			multiplied(cId, innerlupe) = current.getEntry(cId, innerlupe);
		}
	}

	for (cId = 0; cId < getN(); cId++)
	{
		for (innerlupe = 0; innerlupe < Rank.COLNUMBER; innerlupe++)
		{
			fn1 = (*ONEFLOW::ug.c2f)[cId].size();
			n1 = Rank.TempIA[cId];
			for (iFace = 0; iFace < fn1; iFace++)
			{
				//fId = (*ONEFLOW::ug.c2f)[cId][iFace];
				//lc = (*ONEFLOW::ug.lcf)[fId];                                                                     // 面左侧单元
				//rc = (*ONEFLOW::ug.rcf)[fId];                                                                     // 面右侧单元
				offdiag = Rank.TempJA[n1 + iFace];
				if (offdiag > cId && fn2 > 1.0E-16)
				{
					fn2 = (*ONEFLOW::ug.c2f)[offdiag].size();
					n2 = Rank.TempIA[offdiag];
					for (iFace = 0; iFace < fn2; iFace++)
					{
						offdiag2 = Rank.TempJA[n2 + iFace];
						if (offdiag2 == cId)
						{
							multiplied(offdiag, innerlupe) += Rank.TempA[n2 + iFace] * multiplied(offdiag2, innerlupe) / Rank.TempA[n2 + fn2];
						}
					}
				}
			}
		}
	}
	return(multiplied);
}


Solution Preconditioner::solve2(const Solution& current)
{
	Solution multiplied(current);
	int cId;
	int innerlupe;
	for (cId = 0; cId < getN(); cId++)
	{
		for (innerlupe = 0; innerlupe < Rank.COLNUMBER; innerlupe++)
		{
			multiplied(cId, innerlupe) = current.getEntry(cId, innerlupe);
		}
	}
	return(multiplied);
}

#include "util.h"
#include <cmath>
#include <vector>
#include "stdlib.h"
#include "solution.h"
#include "systemSolver.h"
#include<crtdbg.h>

using namespace std;
/** ************************************************************************
 * Update the current approximation to the solution to the linear
 * system. This assumes that the update is created using a GMRES
 * routine, and the upper Hessenberg matrix has been transformed to an
 * upper diagonal matrix already. Note that it changes the values of
 * the values in the coefficients vector, s, which means that the s
 * vector cannot be reused after this without being re-initialized.
 *
 * @return N/A
 ************************************************************************ */
template <class Approximation, class Double >
void Update
(Double** H,         //<! The upper diagonal matrix constructed in the GMRES routine.
	Approximation* x,   //<! The current approximation to the linear system.
	Double** s,          //<! The vector e_1 that has been multiplied by the Givens rotations.
	std::vector<Approximation>* v,  //<! The orthogonal basis vectors for the Krylov subspace.
	int dimension)      //<! The number of vectors in the basis for the Krylov subspace.
{
	// Solve for the coefficients, i.e. solve for c in
	// H*c=s, but we do it in place.
	// 这里分列进行计算
	int lupe;
	int i;
	for (i = 0; i < Rank.COLNUMBER; i++)
	{
		for (lupe = Rank.COLNUMBER * dimension + Rank.COLNUMBER - 1; lupe >= 0; lupe--)
		{
			s[lupe][i] = s[lupe][i] / H[lupe][lupe];
			for (int innerLupe = lupe - 1; innerLupe >= 0; innerLupe--)
			{
				// Subtract off the parts from the upper diagonal of the
				// matrix.
				s[innerLupe][i] -= s[lupe][i] * H[innerLupe][lupe];
			}
		}
	}
	// Finally update the approximation.
	typename std::vector<Approximation>::iterator ptr = v->begin();
	for (int xcol = 0; xcol < Rank.COLNUMBER; xcol++)
	{
		for (int lupe = Rank.RANKNUMBER - 1; lupe >= 0; lupe--)
		{
			for (int m = dimension; m >= 0; m--)
			{
				for (int i = Rank.COLNUMBER - 1; i >= 0; i--)
				{
					(*x)(lupe, xcol) += (*(ptr + m))(lupe, i) * s[Rank.COLNUMBER * m + i][xcol];
				}
			}
		}
	}
}


/** ************************************************************************
 * Implementation of the restarted GMRES algorithm. Follows the
 * algorithm given in the book Templates for the Solution of Linear
 * Systems: Building Blocks for Iterative Methods, 2nd Edition.
 *
 * @return The number of iterations required. Returns zero if it did
 *         not converge.
 ************************************************************************ */
template<class Operation, class Approximation, class Preconditioner, class Double>
int GMRES
(Operation* linearization, //!< Performs the linearization of the PDE on the approximation.
	Approximation* solution,  //!< The approximation to the linear system. (and initial estimate!)
	Approximation* rhs,       //!< the right hand side of the equation to solve.
	Approximation* residual,
	Preconditioner* precond,  //!< The preconditioner used for the linear system.
	int krylovDimension,      //!< The number of vectors to generate in the Krylov subspace.
	int numberRestarts,       //!< Number of times to repeat the GMRES iterations.
	Double tolerance          //!< How small the residual should be to terminate the GMRES iterations.
)
{

	// Allocate the space for the givens rotations, and the upper
	// Hessenburg matrix.
	Double** H = ArrayUtils<Double>::twotensor((krylovDimension + 1) * Rank.COLNUMBER, krylovDimension * Rank.COLNUMBER);
	/*cout << H[1][9] << endl;*/
	// The Givens rotations include the sine and cosine term. The
	// cosine term is in column zero, and the sine term is in column
	// one.
	Double** givens = ArrayUtils<Double>::twotensor(Rank.COLNUMBER * (krylovDimension + 1), krylovDimension * Rank.COLNUMBER);

	Double** s = ArrayUtils<Double>::twotensor((krylovDimension + 1) * Rank.COLNUMBER, Rank.COLNUMBER);
	Double** S = ArrayUtils<Double>::twotensor(Rank.COLNUMBER, Rank.COLNUMBER);  //用来取s矩阵的最后三行
	Double** R = ArrayUtils<Double>::twotensor(Rank.COLNUMBER, Rank.COLNUMBER);
	// Determine the residual and allocate the space for the Krylov
	// subspace.
	std::vector<Approximation> V(krylovDimension + 1,
		Approximation(solution->getN()));
	(*residual) = precond->solve2((*rhs) - (*linearization) * (*solution));

	//残差的第一个列向量
	Double normr1;
	double temp = 0.0;
	for (int ia = 0; ia < Rank.RANKNUMBER; ia++)
	{
		temp += (*residual)(ia, 0) * (*residual)(ia, 0);
	}
	normr1 = sqrt(temp);
	Double normRHS = rhs->norm();
	Double rho = residual->norm();
	Rank.residual = rho;

	// variable for keeping track of how many restarts had to be used.
	int totalRestarts = 0;

	if (normRHS < 1.0E-5)
		normRHS = 1.0;

	// Go through the requisite number of restarts.
	int iteration = 1;
	while ((numberRestarts-- >= 0) && (rho > tolerance * normRHS))
	{

		// The first vector in the Krylov subspace is the normalized residual.
		// The first QR decomposition, and get the vector V[0], V[1], V[2]
		for (int ib = 0; ib < Rank.RANKNUMBER; ib++)
		{
			(V[0])(ib, 0) = (*residual)(ib, 0) * (1.0 / normr1);    //单位化了第一个V向量的第一列，即v1
		}

		R[0][0] = normr1;
		for (int ic = 1; ic < Rank.COLNUMBER; ic++)
		{
			for (int j = 0; j < ic; j++)
			{
				R[j][ic] = 0.0;
				for (int k = 0; k < Rank.RANKNUMBER; k++)
				{
					R[j][ic] += (V[0])(k, j) * (*residual)(k, ic);
				}
				for (int m = 0; m < Rank.RANKNUMBER; m++)
				{
					(V[0])(m, ic) += (V[0])(m, j) * R[j][ic];
				}
			}
			for (int n = 0; n < Rank.RANKNUMBER; n++)
			{
				(V[0])(n, ic) = (*residual)(n, ic) - (V[0])(n, ic);
				R[ic][ic] += V[0](n, ic) * V[0](n, ic);
			}
			R[ic][ic] = sqrt(R[ic][ic]);
			for (int n = 0; n < Rank.RANKNUMBER; n++)
			{
				V[0](n, ic) = V[0](n, ic) * (1 / R[ic][ic]);
			}
		}



		// Need to zero out the s vector in case of restarts
		// initialize the s vector used to estimate the residual.
		for (int lupe = 0; lupe <= krylovDimension * Rank.COLNUMBER + Rank.COLNUMBER - 1; lupe++)
			for (int innerlupe = 0; innerlupe < Rank.COLNUMBER; innerlupe++)
			{
				if (lupe > innerlupe)
					s[lupe][innerlupe] = 0.0;
				else
					s[lupe][innerlupe] = R[lupe][innerlupe];
			}


		// Go through and generate the pre-determined number of vectors
		// for the Krylov subspace.
		for (iteration = 0; iteration < krylovDimension; ++iteration)
		{
			// Get the next entry in the vectors that form the basis for
			// the Krylov subspace.
			V[iteration + 1] = precond->solve2((*linearization) * V[iteration]);
			// Perform the modified Gram-Schmidt method to orthogonalize
			// the new vector.
			int row;
			typename std::vector<Approximation>::iterator ptr = V.begin();
			for (int id = 0; id < Rank.COLNUMBER; id++)
			{
				for (row = 0; row <= Rank.COLNUMBER * iteration + Rank.COLNUMBER - 1; row++)
				{
					for (int n = 0; n < Rank.RANKNUMBER; n++)
					{
						int mod = 0;
						int quo = 0;
						mod = row % Rank.COLNUMBER;
						quo = row / Rank.COLNUMBER;
						H[row][Rank.COLNUMBER * iteration + id] += (V[iteration + 1])(n, id) * (*(ptr + quo))(n, mod);
					}
				}
				//subtract H[row][iteration]*V[row] from the current vector
				for (int row = 0; row <= iteration; row++)
				{
					for (int n = 0; n < Rank.RANKNUMBER; n++)
					{
						double tem = 0.0;
						for (int j = 0; j < Rank.COLNUMBER; j++)
						{
							tem += (*(ptr + row))(n, j) * H[Rank.COLNUMBER * row + j][Rank.COLNUMBER * iteration + id];
						}
						(V[iteration + 1])(n, id) = (V[iteration + 1])(n, id) - tem;
					}
				}
			}

			// QR decomposition to get H[iteration+1][iteration], V[iteration+1]的单位化过程
			std::vector<Approximation> TV(1, Approximation(solution->getN()));   //临时变量
			TV[0] = V[iteration + 1];
			Double normr2;
			double t = 0.0;
			for (int lupe = 0; lupe < Rank.RANKNUMBER; lupe++)
			{
				t += (V[iteration + 1])(lupe, 0) * (V[iteration + 1])(lupe, 0);
			}
			normr2 = sqrt(t);
			for (row = 0; row < Rank.RANKNUMBER; row++)
			{
				(V[iteration + 1])(row, 0) = (V[iteration + 1])(row, 0) * (1.0 / normr2);
			}
			H[Rank.COLNUMBER * (iteration + 1)][Rank.COLNUMBER * iteration] = normr2;
			for (int p = 1; p < Rank.COLNUMBER; p++)
			{
				for (int q = 0; q < p; q++)
				{
					for (int k = 0; k < Rank.RANKNUMBER; k++)
					{
						H[Rank.COLNUMBER * (iteration + 1) + q][Rank.COLNUMBER * iteration + p] += (V[iteration + 1])(k, q) * (TV[0])(k, p);
					}
					for (int m = 0; m < Rank.RANKNUMBER; m++)
					{
						(V[iteration + 1])(m, p) = 0.0;
						(V[iteration + 1])(m, p) += (V[iteration + 1])(m, q) * H[Rank.COLNUMBER * (iteration + 1) + q][Rank.COLNUMBER * iteration + p];
					}
				}
				for (int n = 0; n < Rank.RANKNUMBER; n++)
				{
					(V[iteration + 1])(n, p) = (TV[0])(n, p) - (V[iteration + 1])(n, p);
					H[Rank.COLNUMBER * (iteration + 1) + p][Rank.COLNUMBER * iteration + p] += V[iteration + 1](n, p) * V[iteration + 1](n, p);
				}
				H[Rank.COLNUMBER * (iteration + 1) + p][Rank.COLNUMBER * iteration + p] = sqrt(H[Rank.COLNUMBER * (iteration + 1) + p][Rank.COLNUMBER * iteration + p]);
				for (int n = 0; n < Rank.RANKNUMBER; n++)
				{
					(V[iteration + 1])(n, p) = (V[iteration + 1])(n, p) * (1 / H[Rank.COLNUMBER * (iteration + 1) + p][Rank.COLNUMBER * iteration + p]);
				}
			}
			std::vector<Approximation>(TV).swap(TV);

			// Apply the Givens Rotations to insure that H is
			// an upper diagonal matrix. First apply previous
			// rotations to the current matrix.
			// the first Givens Rotations
			double tmp = 0.0;
			int col;
			for (int lupe = 0; lupe < Rank.COLNUMBER; lupe++)
			{
				for (col = Rank.COLNUMBER * iteration; col < Rank.COLNUMBER * (iteration + 1); col++)
				{
					for (row = 0; row < col; row++)
					{
						int tp = row + Rank.COLNUMBER - 1 - lupe;
						int mp = Rank.COLNUMBER - lupe - 1;
						tmp = givens[tp][2 * (tp - mp)] * H[tp][col]
							+ givens[tp][2 * (tp - mp) + 1] * H[tp + 1][col];
						H[tp + 1][col] = -givens[tp][2 * (tp - mp) + 1] * H[tp][col]
							+ givens[tp][2 * (tp - mp)] * H[tp + 1][col];
						H[tp][col] = tmp;
					}
					// Figure out the next Givens rotation.
					if (H[col + Rank.COLNUMBER - lupe][col] == 0.0)
					{
						// It is already lower diagonal. Just leave it be....
						givens[col + Rank.COLNUMBER - 1 - lupe][col * 2] = 1.0;
						givens[col + Rank.COLNUMBER - 1 - lupe][1 + col * 2] = 0.0;
					}
					else if (fabs(H[col + Rank.COLNUMBER - lupe][col]) > fabs(H[col + Rank.COLNUMBER - 1 - lupe][col]))
					{
						// The off diagonal entry has a larger
						// magnitude. Use the ratio of the
						// diagonal entry over the off diagonal.
						tmp = H[col + Rank.COLNUMBER - 1 - lupe][col] / H[col + Rank.COLNUMBER - lupe][col];
						givens[col + Rank.COLNUMBER - 1 - lupe][1 + col * 2] = 1.0 / sqrt(1.0 + tmp * tmp);
						givens[col + Rank.COLNUMBER - 1 - lupe][col * 2] = tmp * givens[col + Rank.COLNUMBER - 1 - lupe][1 + col * 2];
					}
					else
					{
						// The off diagonal entry has a smaller
						// magnitude. Use the ratio of the off
						// diagonal entry to the diagonal entry.
						tmp = H[col + Rank.COLNUMBER - lupe][col] / H[col + Rank.COLNUMBER - 1 - lupe][col];
						givens[col + Rank.COLNUMBER - 1 - lupe][col * 2] = 1.0 / sqrt(1.0 + tmp * tmp);
						givens[col + Rank.COLNUMBER - 1 - lupe][1 + col * 2] = tmp * givens[col + Rank.COLNUMBER - 1 - lupe][col * 2];
					}
					// Apply the new Givens rotation on the
					// new entry in the uppper Hessenberg matrix.
					tmp = givens[col + Rank.COLNUMBER - 1 - lupe][col * 2] * H[col + Rank.COLNUMBER - 1 - lupe][col] +
						givens[col + Rank.COLNUMBER - 1 - lupe][1 + col * 2] * H[col + Rank.COLNUMBER - lupe][col];
					H[col + Rank.COLNUMBER - lupe][col] = -givens[col + Rank.COLNUMBER - 1 - lupe][1 + col * 2] * H[col + Rank.COLNUMBER - 1 - lupe][col] +
						givens[col + Rank.COLNUMBER - 1 - lupe][col * 2] * H[col + Rank.COLNUMBER - lupe][col];
					H[col + Rank.COLNUMBER - 1 - lupe][col] = tmp;
					if (H[col + Rank.COLNUMBER - lupe][col] < 1e-10)
						H[col + Rank.COLNUMBER - lupe][col] = 0;

					// Finally apply the new Givens rotation on the s
					// vector
					for (int z = 0; z < Rank.COLNUMBER; z++)
					{
						tmp = givens[col + Rank.COLNUMBER - 1 - lupe][col * 2] * s[col + Rank.COLNUMBER - 1 - lupe][z] + givens[col + Rank.COLNUMBER - 1 - lupe][1 + col * 2] * s[col + Rank.COLNUMBER - lupe][z];
						s[col + Rank.COLNUMBER - lupe][z] = -givens[col + Rank.COLNUMBER - 1 - lupe][1 + col * 2] * s[col + Rank.COLNUMBER - 1 - lupe][z] + givens[col + Rank.COLNUMBER - 1 - lupe][col * 2] * s[col + Rank.COLNUMBER - lupe][z];
						s[col + Rank.COLNUMBER - 1 - lupe][z] = tmp;
					}
				}
			}
			//将s向量组的后三行存入到S向量组中，以便后续求其范数
			int lupe, innerlupe;
			for (lupe = 0; lupe < Rank.COLNUMBER; lupe++)
			{
				for (innerlupe = 0; innerlupe < Rank.COLNUMBER; innerlupe++)
				{
					S[lupe][innerlupe] = s[Rank.COLNUMBER * (iteration + 1) + lupe][innerlupe];
				}
			}

			/*int row;*/
			int im;
			int in;
			rho = 0.0;
			for (im = Rank.COLNUMBER - 1; im >= 0; im--)
				for (in = Rank.COLNUMBER - 1; in >= 0; in--)
				{
					rho += S[im][in] * S[im][in];
				}
			rho = sqrt(rho);

			//cout << "iteration:" << iteration << "residual:" << rho << endl;
			/*ofstream file2("residual.txt", ios::app);
			file2 << "residual:" << rho << endl;
			file2.close();*/

			//cout << tolerance * normRHS << endl;
			if (rho < tolerance * normRHS)   //该判别方式在BGMRES中是否有效
			{
				// We are close enough! Update the approximation.
				Update(H, solution, s, &V, iteration);
				(*residual) = precond->solve2((*linearization) * (*solution) - (*rhs));      
				Rank.residual = residual->norm();
				ArrayUtils<double>::deltwotensor(givens);
				ArrayUtils<double>::deltwotensor(H);
				ArrayUtils<double>::deltwotensor(s);
				ArrayUtils<double>::deltwotensor(S);
				ArrayUtils<double>::deltwotensor(R);
				std::vector<Approximation>().swap(V);
				//delete [] V;
				//tolerance = rho/normRHS;
				return(iteration + totalRestarts * krylovDimension);
			}

		} // for(iteration)

		// We have exceeded the number of iterations. Update the
		// approximation and start over.
		totalRestarts += 1;
		Update(H, solution, s, &V, iteration - 1);
		(*residual) = precond->solve2((*linearization) * (*solution) - (*rhs));      //第三个precond
		Rank.residual = residual->norm();

		//cout << "totalRestarts:" << totalRestarts << endl;
	} // while(numberRestarts,rho)


	ArrayUtils<double>::deltwotensor(givens);
	ArrayUtils<double>::deltwotensor(H);
	ArrayUtils<double>::deltwotensor(s);
	ArrayUtils<double>::deltwotensor(S);
	ArrayUtils<double>::deltwotensor(R);
	std::vector<Approximation>().swap(V);
	//delete [] V;
	//tolerance = rho/normRHS;

	if (rho < tolerance * normRHS)
		return(0);

	return(0);
}

/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "UnsteadyConvergence.h"
#include "Ctrl.h"
#include "Iteration.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

void UnsteadyConvergence::Init( int nEqu )
{
    normList.resize( nEqu );
}

void UnsteadyConvergence::Reset()
{
    sum1      = zero;
    sum2      = zero;
    norm0     = zero;
    totalNorm = zero;
    normList  = zero;
}

void UnsteadyConvergence::Accumulate( const RealField & res,
                                      const RealField & q1,
                                      const RealField & q2 )
{
    for ( std::size_t iEqu = 0; iEqu < res.size(); ++ iEqu )
    {
        Real dq_p = res[ iEqu ];
        Real dq_n = q1[ iEqu ] - q2[ iEqu ];

        sum1             += SQR( dq_p );
        sum2             += SQR( dq_n );
        normList[ iEqu ] += SQR( dq_p );
        totalNorm        += SQR( dq_p );
    }
}

void UnsteadyConvergence::Calculate()
{
    if ( ctrl.iConv == 0 )
    {
        conv = sqrt( ABS( sum1 / ( sum2 + SMALL ) ) );
    }
    else if ( ctrl.iConv == 1 )
    {
        if ( Iteration::innerSteps == 1 )
        {
            norm0 = normList[ 0 ];
        }

        conv = normList[ 0 ] / norm0;
    }
    else if ( ctrl.iConv == 2 )
    {
        if ( Iteration::innerSteps == 1 )
        {
            norm0 = totalNorm;
        }

        conv = totalNorm / norm0;
    }
}

EndNameSpace

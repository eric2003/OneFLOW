/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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

#include "UsdData.h"
#include "Ctrl.h"
#include "Iteration.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

UsdBasic usd;

UsdData::UsdData()
{
    ;
}

UsdData::~UsdData()
{
    ;
}

void UsdData::Init()
{
    int nEqu = 1;
    this->InitSub( nEqu );
}

void UsdData::InitSub( int nEqu )
{
    this->InitBasic();
    this->nEqu = nEqu;
    res.resize( nEqu );
    res1.resize( nEqu );
    res2.resize( nEqu );

    prim.resize( nEqu );
    prim1.resize( nEqu );
    prim2.resize( nEqu );

    q.resize( nEqu );
    q1.resize( nEqu );
    q2.resize( nEqu );

    dualtimeRes.resize( nEqu );
    dualtimeSrc.resize( nEqu );

    normList.resize( nEqu );
}

void UsdData::CalcCellDualTimeResidual()
{
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        dualtimeRes[ iEqu ] = resc1 * res [ iEqu ] + 
                               resc2 * res1[ iEqu ] + 
                              resc3 * res2[ iEqu ];
    }
}

void UsdData::CalcCellDualTimeSrc()
{
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        Real dualSrc0 = sc1 * vol  * q [ iEqu ];
        Real dualSrc1 = sc2 * vol1 * q1[ iEqu ];
        Real dualSrc2 = sc3 * vol2 * q2[ iEqu ];

        Real dualSrc = dualSrc0 + dualSrc1 + dualSrc2;

        dualtimeSrc[ iEqu ] = dualSrc;
    }
}

void UsdData::ZeroData()
{
    sum1      = zero;
    sum2      = zero;
    norm0     = zero;
    totalNorm = zero;
    normList  = zero;
}

void UsdData::CalcCellUnsteadyCri()
{
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        Real dq_p =  res[ iEqu ];             // qn+1, p+1 - qn+1, p
        Real dq_n =  q1[ iEqu ] - q2[ iEqu ]; // qn+1, p+1 - qn
        sum1             += SQR( dq_p );
        sum2             += SQR( dq_n );
        normList[ iEqu ] += SQR( dq_p );
        totalNorm        += SQR( dq_p );
    }
}

void UsdData::CalcCvg()
{
    if ( ctrl.iConv == 0 )
    {
        this->conv = sqrt( ABS( sum1 / ( sum2 + SMALL ) ) );
    }
    else if ( ctrl.iConv == 1 )
    {
        if ( Iteration::innerSteps == 1 )
        {
            this->norm0 = this->normList[ 0 ];
        }
        
        this->conv = this->normList[ 0 ] / this->norm0;
    }
    else if ( ctrl.iConv == 2 )
    {
        if ( Iteration::innerSteps == 1 )
        {
            this->norm0 = this->totalNorm;
        }
        this->conv = totalNorm / this->norm0;
    }
}


EndNameSpace
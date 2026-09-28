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
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "UsdData.h"
#include "TimeIntegration.h"

BeginNameSpace( ONEFLOW )



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
    timeIntegration.Init();
    this->nEqu = nEqu;
    res.resize( nEqu );
    res1.resize( nEqu );
    res2.resize( nEqu );

    q.resize( nEqu );
    q1.resize( nEqu );
    q2.resize( nEqu );

    dualtimeRes.resize( nEqu );
    dualtimeSrc.resize( nEqu );

    convergence.Init( nEqu );
}

void UsdData::CalcCellDualTimeResidual()
{
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        dualtimeRes[ iEqu ] = timeIntegration.resc1 * res [ iEqu ] + 
                               timeIntegration.resc2 * res1[ iEqu ] + 
                              timeIntegration.resc3 * res2[ iEqu ];
    }
}

void UsdData::CalcCellDualTimeSrc( Real vol, Real vol1, Real vol2 )
{
    for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
    {
        Real dualSrc0 = timeIntegration.sc1 * vol  * q [ iEqu ];
        Real dualSrc1 = timeIntegration.sc2 * vol1 * q1[ iEqu ];
        Real dualSrc2 = timeIntegration.sc3 * vol2 * q2[ iEqu ];

        Real dualSrc = dualSrc0 + dualSrc1 + dualSrc2;

        dualtimeSrc[ iEqu ] = dualSrc;
    }
}


EndNameSpace

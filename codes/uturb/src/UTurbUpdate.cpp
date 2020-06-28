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

#include "UTurbUpdate.h"
#include "TurbCom.h"
#include "UTurbCom.h"
#include "TurbRhs.h"
#include "UCom.h"
#include "NsCom.h"
#include "UNsCom.h"
#include "NsIdx.h"
#include "Zone.h"
#include "ZoneState.h"
#include "HXMath.h"
#include "Iteration.h"
#include <iostream>
#include <iomanip>
using namespace std;

BeginNameSpace( ONEFLOW )

UTurbUpdate::UTurbUpdate()
{
}

UTurbUpdate::~UTurbUpdate()
{
}

void UTurbUpdate::UpdateFlowField( int sTid )
{
    GetUpdateField( sTid, this->q, this->dq );

    ug.Init();
    uturbf.Init();

    if ( turbcom.nEqu == 1 )
    {
        this->UpdateFlowField1Equ();
    }
    else if ( turbcom.nEqu >= 2 )
    {
        this->UpdateFlowField2Equ();
    }
}

void UTurbUpdate::UpdateFlowField1Equ()
{
    this->UpdateFlowField1EquStd();
}

void UTurbUpdate::UpdateFlowField2Equ()
{
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId = cId;

        this->PrepareData2Equ();

        this->CalcFlowField2Equ();

        this->UpdateValue();
    }
}

void UTurbUpdate::UpdateFlowField1EquStd()
{
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId = cId;

        this->PrepareData1Equ();

        this->CalcFlowField1Equ();

        this->UpdateValue();
    }

}

void UTurbUpdate::UpdateFlowField2EquStd()
{
    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId = cId;

        this->PrepareData2Equ();

        this->CalcFlowField2Equ();

        this->UpdateValue();
    }
}

void UTurbUpdate::UpdateValue()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        ( * uturbf.q )[ iEqu ][ ug.cId  ] = turbcom.q [ iEqu ];
    }
}

void UTurbUpdate::PrepareData1Equ()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbcom.q [ iEqu ] = ( * uturbf.q )[ iEqu ][ ug.cId ];
        turbcom.q0[ iEqu ] = ( * uturbf.q )[ iEqu ][ ug.cId ];
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbcom.dq[ iEqu ] = ( * uturbf.dq )[ iEqu ][ ug.cId ];
    }

    turbcom.rho = ABS( ( * uturbf.q_ns )[ IDX::IR ][ ug.cId ] ) + SMALL;
    turbcom.visl = ( * uturbf.visl )[ 0 ][ ug.cId ];
}

void UTurbUpdate::PrepareData2Equ()
{
    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbcom.q [ iEqu ] = ( * uturbf.q )[ iEqu ][ ug.cId ];
        turbcom.q0[ iEqu ] = ( * uturbf.q )[ iEqu ][ ug.cId ];
    }

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbcom.dq[ iEqu ] = ( * uturbf.dq )[ iEqu ][ ug.cId ];
    }

    turbcom.rho = ABS( ( * uturbf.q_ns )[ IDX::IR ][ ug.cId ] ) + SMALL;
    turbcom.visl = ( * uturbf.visl )[ 0 ][ ug.cId ];
}

void UTurbUpdate::DumpProbeInfo()
{
    cout << setprecision( 3 );
    cout << "Warning : p = " << nscom.prim[ IDX::IP ] << ", r = " << nscom.prim[ IDX::IR ];
    cout << " <-> zid = " << ZoneState::zid << ", cid = " << ug.cId << endl;
}

void UTurbUpdate::SmoothTurbulentPoint()
{
    turbcom.q = 0;

    Real sumV = 0.0;

    int fn = ( * ug.c2f )[ ug.cId ].size();

    for ( int iFace = 0; iFace < fn; ++ iFace )
    {
        int fId = ( * ug.c2f )[ ug.cId ][ iFace ];

        ug.lc = ( * ug.lcf )[ fId ];
        ug.rc = ( * ug.rcf )[ fId ];

        int iNei = ug.lc;
        if ( ug.cId == ug.lc  ) iNei = ug.rc;

        Real volN = one;

        sumV += volN;

        for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
        {
            Real f = ( * unsf.q )[ iEqu ][ iNei ];
            turbcom.q[ iEqu ] += f * volN;
        }
    }

    Real rVol = 1.0 / sumV;

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbcom.q[ iEqu ] *= rVol;
    }
}

EndNameSpace
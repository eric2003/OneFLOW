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

#include "UNsUpdate.h"
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

UNsUpdate::UNsUpdate()
{
}

UNsUpdate::~UNsUpdate()
{
}

void UNsUpdate::UpdateFlowField( int sTid )
{
    GetUpdateField( sTid, this->q, this->dq );

    ug.Init();
    unsf.Init();

    for ( int cId = 0; cId < ug.nCell; ++ cId )
    {
        ug.cId = cId;

        if ( ug.cId == 9 )
        {
            int kkk = 1;
        }

        this->PrepareData();

        this->CalcFlowField();
        //this->CalcFlowFieldHyperSonic();
        //this->CalcFlowFieldHyperSonic_Temperature();

        this->UpdateFlowFieldValue();
    }
}

void UNsUpdate::PrepareData()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.prim [ iEqu ] = ( * unsf.q )[ iEqu ][ ug.cId ];
        nscom.prim0[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.cId ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
    {
        nscom.t [ iEqu ] = ( * unsf.tempr )[ iEqu ][ ug.cId ];
        nscom.t0[ iEqu ] = ( * unsf.tempr )[ iEqu ][ ug.cId ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.dq[ iEqu ] = ( * unsf.dq )[ iEqu ][ ug.cId ];
    }

    nscom.gama = ( * unsf.gama )[ 0 ][ ug.cId ];
}

void UNsUpdate::DumpProbeInfo()
{
    cout << setprecision( 3 );
    cout << "Warning : p = " << nscom.prim[ IDX::IP ] << ", r = " << nscom.prim[ IDX::IR ];
    cout << " <-> zid = " << ZoneState::zid << ", cid = " << ug.cId << endl;
}

void UNsUpdate::SolutionFix()
{
    nscom.prim = 0;

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

        for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
        {
            Real f = ( * unsf.q )[ iEqu ][ iNei ];
            if ( iEqu == IDX::IR || iEqu == IDX::IP )
            {
                f = ABS( f );
            }
            nscom.prim[ iEqu ] += f * volN;
        }
    }

    Real rVol = 1.0 / sumV;

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.prim[ iEqu ] *= rVol;
    }
}

void UNsUpdate::UpdateFlowFieldValue()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        ( * unsf.q )[ iEqu ][ ug.cId ] = nscom.prim[ iEqu ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
    {
        ( * unsf.tempr )[ iEqu ][ ug.cId ] = nscom.t[ iEqu ];
    }
}

EndNameSpace
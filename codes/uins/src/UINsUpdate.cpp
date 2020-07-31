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

#include "UINsUpdate.h"
#include "UCom.h"
#include "INsCom.h"
#include "UINsCom.h"
#include "INsIdx.h"
#include "Zone.h"
#include "ZoneState.h"
#include "HXMath.h"
#include "Iteration.h"
#include <iostream>
#include <iomanip>
using namespace std;

BeginNameSpace( ONEFLOW )

UINsUpdate::UINsUpdate()
{
}

UINsUpdate::~UINsUpdate()
{
}

void UINsUpdate::UpdateFlowField( int sTid )
{
    GetUpdateField( sTid, this->q, this->dq );

    ug.Init();
    uinsf.Init();

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

void UINsUpdate::PrepareData()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.prim [ iEqu ] = ( * uinsf.q )[ iEqu ][ ug.cId ];
        nscom.prim0[ iEqu ] = ( * uinsf.q )[ iEqu ][ ug.cId ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
    {
        nscom.t [ iEqu ] = ( * uinsf.tempr )[ iEqu ][ ug.cId ];
        nscom.t0[ iEqu ] = ( * uinsf.tempr )[ iEqu ][ ug.cId ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.dq[ iEqu ] = ( * uinsf.dq )[ iEqu ][ ug.cId ];
    }

    nscom.gama = ( * uinsf.gama )[ 0 ][ ug.cId ];
}

void UINsUpdate::DumpProbeInfo()
{
    cout << setprecision( 3 );
    cout << "Warning : p = " << nscom.prim[ IIDX::IIP ] << ", r = " << nscom.prim[ IIDX::IIR ];
    cout << " <-> zid = " << ZoneState::zid << ", cid = " << ug.cId << endl;
}

void UINsUpdate::SolutionFix()
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
            Real f = ( * uinsf.q )[ iEqu ][ iNei ];
            if ( iEqu == IIDX::IIR || iEqu == IIDX::IIP )
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

void UINsUpdate::UpdateFlowFieldValue()
{
    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        ( * uinsf.q )[ iEqu ][ ug.cId ] = nscom.prim[ iEqu ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
    {
        ( * uinsf.tempr )[ iEqu ][ ug.cId ] = nscom.t[ iEqu ];
    }
}

EndNameSpace
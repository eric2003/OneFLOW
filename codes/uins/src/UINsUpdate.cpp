/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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
#include "INsInvterm.h"
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

void UINsUpdate::UpdateINsFlowField( int sTid )
{
    GetUpdateField( sTid, this->q, this->dq );

    ug.Init();
    uinsf.Init();

    for ( int cId = 0; cId < ug.nTCell; ++ cId )
    {
        ug.cId = cId;

        if ( ug.cId == 9 )
        {
            int kkk = 1;
        }

        this->PrepareData();

        this->CmpFlowField();
        //this->CmpFlowFieldHyperSonic();
        //this->CmpFlowFieldHyperSonic_Temperature();

        this->UpdateFlowFieldValue();
    }
}

void UINsUpdate::PrepareData()
{
    for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
    {
        inscom.prim [ iEqu ] = ( * uinsf.q )[ iEqu ][ ug.cId ];
        inscom.prim0[ iEqu ] = ( * uinsf.q )[ iEqu ][ ug.cId ];
    }

    for ( int iEqu = 0; iEqu < inscom.nTModel; ++ iEqu )
    {
        inscom.t [ iEqu ] = ( * uinsf.tempr )[ iEqu ][ ug.cId ];
        inscom.t0[ iEqu ] = ( * uinsf.tempr )[ iEqu ][ ug.cId ];
    }

    for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
    {
        inscom.dq[ iEqu ] = ( * uinsf.dq )[ iEqu ][ ug.cId ];
    }

    inscom.gama = ( * uinsf.gama )[ 0 ][ ug.cId ];
}

void UINsUpdate::DumpProbeInfo()
{
    cout << setprecision( 3 );
    cout << "Warning : p = " << inscom.prim[ IIDX::IIP ] << ", r = " << inscom.prim[ IIDX::IIR ];
    cout << " <-> zid = " << ZoneState::zid << ", cid = " << ug.cId << endl;
}

void UINsUpdate::SolutionFix()
{
    inscom.prim = 0;

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

        for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
        {
            Real f = ( * uinsf.q )[ iEqu ][ iNei ];
            if ( iEqu == IIDX::IIR || iEqu == IIDX::IIP )
            {
                f = ABS( f );
            }
            inscom.prim[ iEqu ] += f * volN;
        }
    }

    Real rVol = 1.0 / sumV;

    for ( int iEqu = 0; iEqu < inscom.nTEqu; ++ iEqu )
    {
        inscom.prim[ iEqu ] *= rVol;
    }
}

void UINsUpdate::UpdateFlowFieldValue()
{
	(*uinsf.q)[IIDX::IIU][ug.cId] = iinv.up[ug.cId];
	(*uinsf.q)[IIDX::IIV][ug.cId] = iinv.vp[ug.cId];
	(*uinsf.q)[IIDX::IIW][ug.cId] = iinv.wp[ug.cId];
	(*uinsf.q)[IIDX::IIP][ug.cId] = iinv.pc[ug.cId];

	//for (int iEqu = 0; iEqu < inscom.nTEqu; ++iEqu)
	//{
		//(*uinsf.q)[iEqu][ug.cId] = inscom.prim[iEqu];
	//}


    for ( int iEqu = 0; iEqu < inscom.nTModel; ++ iEqu )
    {
        ( * uinsf.tempr )[ iEqu ][ ug.cId ] = inscom.t[ iEqu ];
    }
}

EndNameSpace
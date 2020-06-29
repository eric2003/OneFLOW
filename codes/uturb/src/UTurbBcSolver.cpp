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

#include "UTurbBcSolver.h"
#include "UTurbCom.h"
#include "TurbCom.h"
#include "UBcSolver.h"
#include "Ctrl.h"
#include "BcData.h"
#include "NsCom.h"
#include "UCom.h"
#include "UNsCom.h"
#include "UnsGrid.h"
#include "DataBase.h"
#include "Zone.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "Boundary.h"
#include "HXStd.h"

BeginNameSpace( ONEFLOW )

UTurbBcSolver::UTurbBcSolver()
{
}

UTurbBcSolver::~UTurbBcSolver()
{
}

void UTurbBcSolver::Init()
{
    ug.Init();
    uturbf.Init();
    turbcom.Init();
}

void UTurbBcSolver::CalcBc()
{
    ug.nRegion = ug.bcRecord->bcInfo->bcType.size();

    for ( int ir = 0; ir < ug.nRegion; ++ ir )
    {
        ug.ir = ir;
        ug.bctype = ug.bcRecord->bcInfo->bcType[ ir ];
        ug.nRBFace = ug.bcRecord->bcInfo->bcFace[ ir ].size();
        this->SetBc();

        this->CalcBcRegion();
    }
}

void UTurbBcSolver::SetId( int bcfId )
{
    ug.bcfId = bcfId;

    BcInfo * bcInfo = ug.bcRecord->bcInfo;

    ug.fId = bcInfo->bcFace[ ug.ir ][ bcfId ];
    ug.bcr = bcInfo->bcRegion[ ug.ir ][ bcfId ];

    ug.bcdtkey = bcInfo->bcdtkey[ ug.ir ][ bcfId ];

    ug.lc = ( * ug.lcf )[ ug.fId ];
    ug.rc = ( * ug.rcf )[ ug.fId ];

    turbcom.bcdtkey = 0;
    if ( ug.bcr == -1 ) return;
    int dd = bcdata.r2d[ ug.bcr ];
    if ( dd != - 1 )
    {
        turbcom.bcdtkey = 1;
        turbcom.bcflow = & bcdata.dataList[ dd ];
    }

}

void UTurbBcSolver::CalcBcRegion()
{
    for ( int ibc = 0; ibc < ug.nRBFace; ++ ibc )
    {
        this->SetId( ibc );

        this->PrepareData();

        this->CalcFaceBc();

        this->UpdateBc();
    }
}

void UTurbBcSolver::UpdateBc()
{
    if ( ! this->updateFlag ) return;

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        ( * uturbf.q )[ iEqu ][ ug.rc ] = turbcom.primt1[ iEqu ];
    }
}

void UTurbBcSolver::PrepareData()
{
    gcom.xfn   = ( * ug.xfn   )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn   )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn   )[ ug.fId ];

    gcom.vfx   = ( * ug.vfx   )[ ug.fId ];
    gcom.vfy   = ( * ug.vfy   )[ ug.fId ];
    gcom.vfz   = ( * ug.vfz   )[ ug.fId ];

    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    gcom.xcc1 = ( * ug.xcc )[ ug.lc ];
    gcom.ycc1 = ( * ug.ycc )[ ug.lc ];
    gcom.zcc1 = ( * ug.zcc )[ ug.lc ];

    gcom.xcc2 = ( * ug.xcc )[ ug.rc ];
    gcom.ycc2 = ( * ug.ycc )[ ug.rc ];
    gcom.zcc2 = ( * ug.zcc )[ ug.rc ];

    gcom.xfc =  ( * ug.xfc )[ ug.fId ];
    gcom.yfc =  ( * ug.yfc )[ ug.fId ];
    gcom.zfc =  ( * ug.zfc )[ ug.fId ];

    turbcom.dist  = ( * uturbf.dist )[ ug.lc ];

    for ( int iEqu = 0; iEqu < turbcom.nEqu; ++ iEqu )
    {
        turbcom.prims1[ iEqu ] = ( * uturbf.q )[ iEqu ][ ug.lc ];
        turbcom.prims2[ iEqu ] = ( * uturbf.q )[ iEqu ][ ug.lc ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        turbcom.ns_prims1[ iEqu ] = ( * uturbf.q_ns )[ iEqu ][ ug.lc ];
        turbcom.ns_prims2[ iEqu ] = ( * uturbf.q_ns )[ iEqu ][ ug.lc ];
    }
}

EndNameSpace
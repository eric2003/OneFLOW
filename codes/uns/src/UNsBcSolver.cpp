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

#include "UNsBcSolver.h"
#include "UBcSolver.h"
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

UNsBcSolver::UNsBcSolver()
{
}

UNsBcSolver::~UNsBcSolver()
{
}

void UNsBcSolver::Init()
{
    ug.Init();
    unsf.Init();
}

void UNsBcSolver::CalcBc()
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

void UNsBcSolver::SetId( int bcfId )
{
    ug.bcfId = bcfId;

    BcInfo * bcInfo = ug.bcRecord->bcInfo;

    ug.fId = bcInfo->bcFace[ ug.ir ][ bcfId ];
    ug.bcr = bcInfo->bcRegion[ ug.ir ][ bcfId ];

    ug.bcdtkey = bcInfo->bcdtkey[ ug.ir ][ bcfId ];

    ug.lc = ( * ug.lcf )[ ug.fId ];
    ug.rc = ( * ug.rcf )[ ug.fId ];

    nscom.bcdtkey = 0;
    if ( ug.bcr == -1 ) return; //interface
    int dd = bcdata.r2d[ ug.bcr ];
    if ( dd != - 1 )
    {
        nscom.bcdtkey = 1;
        nscom.bcflow = & bcdata.dataList[ dd ];
    }

}

void UNsBcSolver::CalcBcRegion()
{
    for ( int ibc = 0; ibc < ug.nRBFace; ++ ibc )
    {
        this->SetId( ibc );

        this->PrepareData();

        this->CalcFaceBc();

        this->UpdateBc();
    }
}

void UNsBcSolver::UpdateBc()
{
    if ( ! this->updateFlag ) return;

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        ( * unsf.q )[ iEqu ][ ug.rc ] = nscom.primt1[ iEqu ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        ( * unsf.bc_q )[ iEqu ][ ug.fId ] = nscom.prim[ iEqu ];
    }
}

void UNsBcSolver::PrepareData()
{
    gcom.xfn   = ( * ug.xfn   )[ ug.fId ];
    gcom.yfn   = ( * ug.yfn   )[ ug.fId ];
    gcom.zfn   = ( * ug.zfn   )[ ug.fId ];

    gcom.vfx   = ( * ug.vfx   )[ ug.fId ];
    gcom.vfy   = ( * ug.vfy   )[ ug.fId ];
    gcom.vfz   = ( * ug.vfz   )[ ug.fId ];

    gcom.vfn   = ( * ug.vfn   )[ ug.fId ];
    gcom.farea = ( * ug.farea )[ ug.fId ];

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.q1[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.lc ];
        nscom.q2[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.lc ];
    }

    nscom.gama1 = ( * unsf.gama )[ 0 ][ ug.lc ];
    nscom.gama2 = ( * unsf.gama )[ 0 ][ ug.lc ];

    gcom.xcc1 = ( * ug.xcc )[ ug.lc ];
    gcom.ycc1 = ( * ug.ycc )[ ug.lc ];
    gcom.zcc1 = ( * ug.zcc )[ ug.lc ];

    gcom.xcc2 = ( * ug.xcc )[ ug.rc ];
    gcom.ycc2 = ( * ug.ycc )[ ug.rc ];
    gcom.zcc2 = ( * ug.zcc )[ ug.rc ];

    gcom.xfc =  ( * ug.xfc )[ ug.fId ];
    gcom.yfc =  ( * ug.yfc )[ ug.fId ];
    gcom.zfc =  ( * ug.zfc )[ ug.fId ];

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.prims1[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.lc ];
        nscom.prims2[ iEqu ] = ( * unsf.q )[ iEqu ][ ug.lc ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
    {
        nscom.ts1[ iEqu ] = ( * unsf.tempr )[ iEqu ][ ug.lc ];
        nscom.ts2[ iEqu ] = ( * unsf.tempr )[ iEqu ][ ug.lc ];
    }
}

EndNameSpace
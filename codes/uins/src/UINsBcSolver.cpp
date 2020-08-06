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

#include "UINsBcSolver.h"
#include "UBcSolver.h"
#include "BcData.h"
#include "INsCom.h"
#include "UCom.h"
#include "UINsCom.h"
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

UINsBcSolver::UINsBcSolver()
{
}

UINsBcSolver::~UINsBcSolver()
{
}

void UINsBcSolver::Init()
{
    ug.Init();
    uinsf.Init();
}

void UINsBcSolver::CalcBc()
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

void UINsBcSolver::SetId( int bcfId )
{
    ug.bcfId = bcfId;

    BcInfo * bcInfo = ug.bcRecord->bcInfo;

    ug.fId = bcInfo->bcFace[ ug.ir ][ bcfId ];
    ug.bcNameId = bcInfo->bcNameId[ ug.ir ][ bcfId ];

    ug.lc = ( * ug.lcf )[ ug.fId ];
    ug.rc = ( * ug.rcf )[ ug.fId ];

    nscom.bcdtkey = 0;
    if ( ug.bcNameId == -1 ) return; //interface
    int dd = ins_bc_data.r2d[ ug.bcNameId ];
    if ( dd != - 1 )
    {
        nscom.bcdtkey = 1;
        nscom.bcflow = & ins_bc_data.dataList[ dd ];
    }

}

void UINsBcSolver::CalcBcRegion()
{
    for ( int ibc = 0; ibc < ug.nRBFace; ++ ibc )
    {
        this->SetId( ibc );

        this->PrepareData();

        this->CalcFaceBc();

        this->UpdateBc();
    }
}

void UINsBcSolver::UpdateBc()
{
    if ( ! this->updateFlag ) return;

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        ( * uinsf.q )[ iEqu ][ ug.rc ] = nscom.primt1[ iEqu ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        ( * uinsf.bc_q )[ iEqu ][ ug.fId ] = nscom.prim[ iEqu ];
    }
}

void UINsBcSolver::PrepareData()
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
        nscom.q1[ iEqu ] = ( * uinsf.q )[ iEqu ][ ug.lc ];
        nscom.q2[ iEqu ] = ( * uinsf.q )[ iEqu ][ ug.lc ];
    }

    nscom.gama1 = ( * uinsf.gama )[ 0 ][ ug.lc ];
    nscom.gama2 = ( * uinsf.gama )[ 0 ][ ug.lc ];

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
        nscom.prims1[ iEqu ] = ( * uinsf.q )[ iEqu ][ ug.lc ];
        nscom.prims2[ iEqu ] = ( * uinsf.q )[ iEqu ][ ug.lc ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
    {
        nscom.ts1[ iEqu ] = ( * uinsf.tempr )[ iEqu ][ ug.lc ];
        nscom.ts2[ iEqu ] = ( * uinsf.tempr )[ iEqu ][ ug.lc ];
    }
}

EndNameSpace
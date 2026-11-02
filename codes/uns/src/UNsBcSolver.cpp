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
#include "StageProfiler.h"

#include <chrono>

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
    this->CacheFieldPointers();
}

void UNsBcSolver::CacheFieldPointers()
{
    qData.resize( unsf.q->GetNEqu() );
    for ( int iEqu = 0; iEqu < static_cast< int >( qData.size() ); ++ iEqu )
    {
        qData[ iEqu ] = ( * unsf.q )[ iEqu ].data();
    }

    bcQData.resize( unsf.bc_q->GetNEqu() );
    for ( int iEqu = 0; iEqu < static_cast< int >( bcQData.size() ); ++ iEqu )
    {
        bcQData[ iEqu ] = ( * unsf.bc_q )[ iEqu ].data();
    }

    temprData.resize( unsf.tempr->GetNEqu() );
    for ( int iEqu = 0; iEqu < static_cast< int >( temprData.size() ); ++ iEqu )
    {
        temprData[ iEqu ] = ( * unsf.tempr )[ iEqu ].data();
    }

    gamaData = ( * unsf.gama )[ 0 ].data();

    leftCellData = ug.lcf->data();
    rightCellData = ug.rcf->data();
    xfnData = ug.xfn->data();
    yfnData = ug.yfn->data();
    zfnData = ug.zfn->data();
    vfxData = ug.vfx->data();
    vfyData = ug.vfy->data();
    vfzData = ug.vfz->data();
    vfnData = ug.vfn->data();
    faceAreaData = ug.farea->data();
    faceCenterXData = ug.xfc->data();
    faceCenterYData = ug.yfc->data();
    faceCenterZData = ug.zfc->data();
    cellCenterXData = ug.xcc->data();
    cellCenterYData = ug.ycc->data();
    cellCenterZData = ug.zcc->data();
}

void UNsBcSolver::CalcBc()
{
    ug.nRegion = ug.bcRecord->GetBcInfo()->bcType.size();

    for ( int ir = 0; ir < ug.nRegion; ++ ir )
    {
        ug.ir = ir;
        ug.bctype = ug.bcRecord->GetBcInfo()->bcType[ ir ];
        ug.nRBFace = ug.bcRecord->GetBcInfo()->bcFace[ ir ].size();
        this->SetBc();

        this->CalcBcRegion();
    }
}

void UNsBcSolver::SetId( int bcfId )
{
    ug.bcfId = bcfId;

    BcInfo * bcInfo = ug.bcRecord->GetBcInfo();

    ug.fId = bcInfo->bcFace[ ug.ir ][ bcfId ];
    ug.bcNameId = bcInfo->bcNameId[ ug.ir ][ bcfId ];

    ug.lc = leftCellData[ ug.fId ];
    ug.rc = rightCellData[ ug.fId ];

    nscom.bcdtkey = 0;
    if ( ug.bcNameId == -1 ) return; //interface
    int dd = ns_bc_data.r2d[ ug.bcNameId ];
    if ( dd != - 1 )
    {
        nscom.bcdtkey = 1;
        nscom.bcflow = & ns_bc_data.dataList[ dd ];
    }

}

void UNsBcSolver::CalcBcRegion()
{
    static const bool detailEnabled = StageProfiler::Enabled();
    double detailMs[ 4 ] = { 0.0, 0.0, 0.0, 0.0 };

    for ( int ibc = 0; ibc < ug.nRBFace; ++ ibc )
    {
        if ( detailEnabled )
        {
            auto start = std::chrono::steady_clock::now();
            this->SetId( ibc );
            detailMs[ 0 ] += std::chrono::duration< double, std::milli >(
                std::chrono::steady_clock::now() - start ).count();

            start = std::chrono::steady_clock::now();
            this->PrepareData();
            detailMs[ 1 ] += std::chrono::duration< double, std::milli >(
                std::chrono::steady_clock::now() - start ).count();

            start = std::chrono::steady_clock::now();
            this->CalcFaceBc();
            detailMs[ 2 ] += std::chrono::duration< double, std::milli >(
                std::chrono::steady_clock::now() - start ).count();

            start = std::chrono::steady_clock::now();
            this->UpdateBc();
            detailMs[ 3 ] += std::chrono::duration< double, std::milli >(
                std::chrono::steady_clock::now() - start ).count();
        }
        else
        {
            this->SetId( ibc );
            this->PrepareData();
            this->CalcFaceBc();
            this->UpdateBc();
        }
    }

    if ( detailEnabled )
    {
        StageProfiler::Record( "boundary_set_id", detailMs[ 0 ] );
        StageProfiler::Record( "boundary_prepare_data", detailMs[ 1 ] );
        StageProfiler::Record( "boundary_calc_face_bc", detailMs[ 2 ] );
        StageProfiler::Record( "boundary_update_bc", detailMs[ 3 ] );
    }
}

void UNsBcSolver::UpdateBc()
{
    if ( ! this->updateFlag ) return;

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        qData[ iEqu ][ ug.rc ] = nscom.primt1[ iEqu ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        bcQData[ iEqu ][ ug.fId ] = nscom.prim[ iEqu ];
    }
}

void UNsBcSolver::PrepareData()
{
    gcom.xfn   = xfnData[ ug.fId ];
    gcom.yfn   = yfnData[ ug.fId ];
    gcom.zfn   = zfnData[ ug.fId ];

    gcom.vfx   = vfxData[ ug.fId ];
    gcom.vfy   = vfyData[ ug.fId ];
    gcom.vfz   = vfzData[ ug.fId ];

    gcom.vfn   = vfnData[ ug.fId ];
    gcom.farea = faceAreaData[ ug.fId ];

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.q1[ iEqu ] = qData[ iEqu ][ ug.lc ];
        nscom.q2[ iEqu ] = qData[ iEqu ][ ug.lc ];
    }

    nscom.gama1 = gamaData[ ug.lc ];
    nscom.gama2 = gamaData[ ug.lc ];
    nscom.gama  = half * ( nscom.gama1 + nscom.gama2 );

    gcom.xcc1 = cellCenterXData[ ug.lc ];
    gcom.ycc1 = cellCenterYData[ ug.lc ];
    gcom.zcc1 = cellCenterZData[ ug.lc ];

    gcom.xcc2 = cellCenterXData[ ug.rc ];
    gcom.ycc2 = cellCenterYData[ ug.rc ];
    gcom.zcc2 = cellCenterZData[ ug.rc ];

    gcom.xfc = faceCenterXData[ ug.fId ];
    gcom.yfc = faceCenterYData[ ug.fId ];
    gcom.zfc = faceCenterZData[ ug.fId ];

    for ( int iEqu = 0; iEqu < nscom.nTEqu; ++ iEqu )
    {
        nscom.prims1[ iEqu ] = qData[ iEqu ][ ug.lc ];
        nscom.prims2[ iEqu ] = qData[ iEqu ][ ug.lc ];
    }

    for ( int iEqu = 0; iEqu < nscom.nTModel; ++ iEqu )
    {
        nscom.ts1[ iEqu ] = temprData[ iEqu ][ ug.lc ];
        nscom.ts2[ iEqu ] = temprData[ iEqu ][ ug.lc ];
    }
}

EndNameSpace

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
#include "NsIdx.h"
#include "HipProductionPolicy.h"
#ifdef ONEFLOW_ENABLE_HIP
#include "HipTurbulenceBackend.h"
#endif

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
    if ( this->CalcBcHip() ) return;

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

bool UTurbBcSolver::CalcBcHip()
{
#ifdef ONEFLOW_ENABLE_HIP
    if ( turbcom.nEqu != 1 ) return false;
    if ( ! HipProductionFeatureEnabled( "ONEFLOW_ENABLE_UTURB_HIP_BOUNDARY" ) )
    {
        return false;
    }
    UnsGrid * grid = Zone::GetUnsGrid();
    if ( grid == nullptr || ug.nBFaces <= 0 ) return false;
    if ( uturbf.q == nullptr || uturbf.meanFlow.q == nullptr
         || ug.bcRecord == nullptr )
    {
        return false;
    }

    TurbSaBoundaryView view;
    view.nBoundaryFaces = ug.nBFaces;
    view.nCells = ug.nCells;
    view.nGhostCells = ug.nTCell - ug.nCells;
    view.leftCell = ug.lcf->data();
    view.rightCell = ug.rcf->data();
    view.boundaryType = ug.bcRecord->bcType.data();
    view.xNormal = ug.xfn->data();
    view.yNormal = ug.yfn->data();
    view.zNormal = ug.zfn->data();
    view.meshVelocityNormal = ug.vfn->data();
    view.faceOuterNormal = nscom.faceOuterNormal;
    view.workingVariable = ( * uturbf.q )[ ISA ].data();
    view.density = ( * uturbf.meanFlow.q )[ IDX::IR ].data();
    view.pressure = ( * uturbf.meanFlow.q )[ IDX::IP ].data();
    view.velocity[ 0 ] = ( * uturbf.meanFlow.q )[ IDX::IU ].data();
    view.velocity[ 1 ] = ( * uturbf.meanFlow.q )[ IDX::IV ].data();
    view.velocity[ 2 ] = ( * uturbf.meanFlow.q )[ IDX::IW ].data();
    view.inflowWorkingVariable = turbcom.inflow[ ISA ];
    view.inflowDensity = nscom.inflow[ IDX::IR ];
    view.inflowPressure = nscom.inflow[ IDX::IP ];
    view.inflowVelocity[ 0 ] = nscom.inflow[ IDX::IU ];
    view.inflowVelocity[ 1 ] = nscom.inflow[ IDX::IV ];
    view.inflowVelocity[ 2 ] = nscom.inflow[ IDX::IW ];
    view.gammaRef = nscom.gama_ref;
    // nscom.gama is only assigned inside a face loop (UNsInvFlux.cpp:1191), so
    // it can still be zero when the boundary stage runs.  The farfield decision
    // needs a sound speed, and gama_ref is the same constant on every path that
    // matters here: the SA boundary itself uses nscom.gama for the interior
    // sound speed, and the host falls back to the reference value whenever the
    // per-face average has not been computed.  Use gama_ref when gama is unset
    // so the device view is valid; the two agree for the constant-gamma gas
    // this solver runs.
    view.gamma = nscom.gama > 0.0 ? nscom.gama : nscom.gama_ref;
    view.cacheKey = grid;

    HipTurbulenceBackend::Shared().CalcBoundary( view );
    return true;
#else
    return false;
#endif
}

void UTurbBcSolver::SetId( int bcfId )
{
    ug.bcfId = bcfId;

    BcInfo * bcInfo = ug.bcRecord->GetBcInfo();

    ug.fId = bcInfo->bcFace[ ug.ir ][ bcfId ];
    ug.bcNameId = bcInfo->bcNameId[ ug.ir ][ bcfId ];


    ug.lc = ( * ug.lcf )[ ug.fId ];
    ug.rc = ( * ug.rcf )[ ug.fId ];

    turbcom.bcdtkey = 0;
    if ( ug.bcNameId == -1 ) return;
    int dd = turb_bc_data.r2d[ ug.bcNameId ];
    if ( dd != - 1 )
    {
        turbcom.bcdtkey = 1;
        turbcom.bcflow = & turb_bc_data.dataList[ dd ];
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
        turbcom.ns_prims1[ iEqu ] = ( * uturbf.meanFlow.q )[ iEqu ][ ug.lc ];
        turbcom.ns_prims2[ iEqu ] = ( * uturbf.meanFlow.q )[ iEqu ][ ug.lc ];
    }
}

EndNameSpace

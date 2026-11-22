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

#include "UTurbGrad.h"
#include "UCom.h"
#include "TurbCom.h"
#include "DataBase.h"
#include "FieldWrap.h"
#include "SolverDef.h"
#include "FieldManager.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "UnsGrid.h"
#include "Zone.h"
#include "FaceMesh.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include "HipProductionPolicy.h"
#ifdef ONEFLOW_ENABLE_HIP
#include "HipTurbulenceBackend.h"
#endif


BeginNameSpace( ONEFLOW )

UTurbGrad uturb_grad;

UTurbGrad::UTurbGrad()
{
    ;
}

UTurbGrad::~UTurbGrad()
{
    ;
}

void UTurbGrad::Init()
{
    UnsGrid * grid = Zone::GetUnsGrid();

    // Pin solver identity (do not rely on ambient SolverState).
    const int sid = TURB_SOLVER;
    name  = FieldHome::GetGradName( sid, 0 );
    namex = FieldHome::GetGradName( sid, 1 );
    namey = FieldHome::GetGradName( sid, 2 );
    namez = FieldHome::GetGradName( sid, 3 );

    q    = FieldHome::GetSolutionField( grid, sid );
    dqdx = FieldHome::GetGradField( grid, sid, 1 );
    dqdy = FieldHome::GetGradField( grid, sid, 2 );
    dqdz = FieldHome::GetGradField( grid, sid, 3 );

    bdqdx = FieldHome::GetBcGradField( grid, sid, 0 );
    bdqdy = FieldHome::GetBcGradField( grid, sid, 1 );
    bdqdz = FieldHome::GetBcGradField( grid, sid, 2 );

    this->nEqu = turbcom.nEqu;

    this->istore = 1;
}

bool UTurbGrad::CalcGradHip()
{
#ifdef ONEFLOW_ENABLE_HIP
    // One-equation closure only: the backend entry carries a single working
    // variable plane, while a two-equation model would have to supply k and
    // omega separately.
    if ( turbcom.nEqu != 1 ) return false;
    if ( ! HipProductionFeatureEnabled( "ONEFLOW_ENABLE_UTURB_HIP_GRADIENT" ) )
    {
        return false;
    }
    if ( q == nullptr || dqdx == nullptr || dqdy == nullptr || dqdz == nullptr )
    {
        return false;
    }
    UnsGrid * grid = Zone::GetUnsGrid();
    if ( grid == nullptr ) return false;
    ug.Init();

    TurbSaGradientView view;
    view.nCells = ug.nCells;
    view.nGhostCells = ug.nTCell - ug.nCells;
    view.nFaces = ug.nFaces;
    view.nBoundaryFaces = ug.nBFaces;
    view.workingVariable = ( * q )[ 0 ].data();
    view.xFace = ug.xfc->data();
    view.yFace = ug.yfc->data();
    view.zFace = ug.zfc->data();
    view.xNormal = ug.xfn->data();
    view.yNormal = ug.yfn->data();
    view.zNormal = ug.zfn->data();
    view.faceArea = ug.farea->data();
    view.xCell = ug.xcc->data();
    view.yCell = ug.ycc->data();
    view.zCell = ug.zcc->data();
    view.cellVolume = ug.cvol->data();
    view.leftCell = ug.lcf->data();
    view.rightCell = ug.rcf->data();
    view.gradientX = ( * dqdx )[ 0 ].data();
    view.gradientY = ( * dqdy )[ 0 ].data();
    view.gradientZ = ( * dqdz )[ 0 ].data();
    view.cacheKey = grid;

    HipTurbulenceBackend::Shared().CalcGradient( view );
    // The host source term and the diffusive flux read these planes on the
    // next call, and CalcGradient already synchronises and copied the device
    // result into them.
    return true;
#else
    return false;
#endif
}

EndNameSpace

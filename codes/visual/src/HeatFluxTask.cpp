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

#include "HeatFluxTask.h"
#include "HeatFlux.h"
#include "FaceJoint.h"
#include "FileUtil.h"
#include "Prj.h"
#include "WallVisual.h"
#include "AeroForceTask.h"
#include "WallDist.h"
#include "Boundary.h"
#include "BcRecord.h"
#include "NsCtrl.h"
#include "PointSearch.h"
#include "Zone.h"
#include "ZoneState.h"
#include "HXMath.h"
#include "Tolerence.h"
#include "HXCgns.h"
#include "Parallel.h"
#include "ActionState.h"
#include "DataBook.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "FaceTopo.h"
#include "NodeMesh.h"
#include <fstream>
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

HeatFluxTask::HeatFluxTask()
{
    wallManager = new FaceJointManager();
}

HeatFluxTask::~HeatFluxTask()
{
    delete wallManager;
}

void HeatFluxTask::Run()
{
    if ( ns_ctrl.isowallbc == 0 ) return;
    this->AllocVariable();
    this->CollectWallFaceNode();
    this->ConstructPointIndex();
    this->CollectWallFaceValue();
    this->CalcNodeValue();
    this->VisualizeWallNodeValue();
}

void HeatFluxTask::AllocVariable()
{
    for ( int zId = 0; zId < ZoneState::nZones; ++ zId )
    {
        ZoneState::zid = zId;
        wallManager->patch.push_back( new FaceJoint() );
    }
}

void HeatFluxTask::CollectWallFaceNode()
{
    ActionState::dataBook = this->dataBook;

    for ( int zId = 0; zId < ZoneState::nZones; ++ zId )
    {
        ZoneState::zid = zId;

        if ( ZoneState::pid[ zId ] == Parallel::pid )
        {
            ONEFLOW::CollectWallFaceNode();
        }

        HXBcast( ActionState::dataBook, ZoneState::pid[ zId ] );

        AddWallFaceNode( wallManager, zId );
    }
}

void HeatFluxTask::CollectWallFaceValue()
{
    ActionState::dataBook = this->dataBook;

    for ( int zId = 0; zId < ZoneState::nZones; ++ zId )
    {
        ZoneState::zid = zId;

        if ( ZoneState::pid[ zId ] == Parallel::pid )
        {
            ONEFLOW::CollectWallFaceValue();
        }

        HXBcast( ActionState::dataBook, ZoneState::pid[ zId ] );

        AddWallFaceValue( wallManager, zId );
    }
}

void HeatFluxTask::ConstructPointIndex()
{
    wallManager->ConstructPointIndex();
}

void HeatFluxTask::CalcNodeValue()
{
    wallManager->CalcNodeValue();
}

void HeatFluxTask::VisualizeWallNodeValue()
{
    fstream file;
    ONEFLOW::OpenPrjFile( file, ctrl.heatfluxFile, ios_base::out );

    int numberOfSubData = wallManager->patch.size();
    int iCount = 0;
    for ( int iData = 0; iData < numberOfSubData; ++ iData )
    {
        FaceJoint * basicWall = wallManager->patch[ iData ];
        if ( basicWall->isValid ) iCount ++;
        basicWall->Visual( file );
    }
    CloseFile( file );
}

void CollectWallFaceNode()
{
    Grid * gridIn = Zone::GetGrid();
    UnsGrid * grid = UnsGridCast( gridIn );

    int nSolidCell = GetNSolidCell( grid );

    ActionState::dataBook->MoveToBegin();

    HXWrite( ActionState::dataBook, nSolidCell );

    if ( nSolidCell <= 0 ) return;

    WallStructure::PointLink ptLink;

    RealField & x = grid->nodeMesh->xN;
    RealField & y = grid->nodeMesh->yN;
    RealField & z = grid->nodeMesh->zN;

    LinkField & f2n = grid->faceTopo->f2n;
    IntField & bcType = grid->faceTopo->bcManager->bcRecord->bcType;

    int nBFace = bcType.size();

    for ( int iFace = 0; iFace < nBFace; ++ iFace )
    {
        int bc_type = bcType[ iFace ];

        if ( bc_type == BC::SOLID_SURFACE )
        {
            WallStructure::PointField ptList;
            int nNode = f2n[ iFace ].size();
            for ( int iNode = 0; iNode < nNode; ++ iNode )
            {
                int index = f2n[ iFace ][ iNode ];
                Real x0 = x[ index ];
                Real y0 = y[ index ];
                Real z0 = z[ index ];

                WallStructure::PointType pt( x0, y0, z0 );
                ptList.push_back( pt );
            }
            ptLink.push_back( ptList );
        }
    }

    HXWrite( ActionState::dataBook, ptLink );
}

void AddWallFaceNode( FaceJointManager * walldata, int iZone )
{
    ActionState::dataBook->MoveToBegin();
    int nSolidCell;
    HXRead( ActionState::dataBook, nSolidCell );
    
    WallStructure::PointLink ptLink;

    ptLink.resize( nSolidCell );
    HXRead( ActionState::dataBook, ptLink );

    if ( nSolidCell > 0 )
    {
        FaceJoint * global = walldata->global;
        FaceJoint * local = walldata->patch[ iZone ];
        local->isValid = true;
        global->isValid = true;
        local->AddFacePoint( nSolidCell, ptLink );
        global->AddFacePoint( nSolidCell, ptLink );
    }
}

void AddWallFaceValue( FaceJointManager * walldata, int iZone )
{
    ActionState::dataBook->MoveToBegin();
    int nSolidCell;
    HXRead( ActionState::dataBook, nSolidCell );

    RealField fcv;

    fcv.resize( nSolidCell );
    HXRead( ActionState::dataBook, fcv );

    if ( nSolidCell > 0 )
    {
        FaceJoint * global = walldata->global;
        FaceJoint * local  = walldata->patch[ iZone ];

        global->AddFaceCenterValue( nSolidCell, fcv );
        local->AddFaceCenterValue( nSolidCell, fcv );
    }
}

void CollectWallFaceValue()
{
    ActionState::dataBook->MoveToBegin();
    ActionState::dataBook->ReSize( 0 );

    Grid * gridIn = Zone::GetGrid();
    UnsGrid * grid = UnsGridCast( gridIn );

    int nSolidCell = GetNSolidCell( grid );

    HXWrite( ActionState::dataBook, nSolidCell );

    if ( nSolidCell == 0 ) return;

    int zId = ZoneState::zid;

    SurfaceValue * heat_sur = heat_flux.heatflux[ zId ];
    SurfaceValue * fric_sur = heat_flux.fricflux[ zId ];

    RealField & hf = * heat_sur->var;

    HXWrite( ActionState::dataBook, hf );
}



EndNameSpace
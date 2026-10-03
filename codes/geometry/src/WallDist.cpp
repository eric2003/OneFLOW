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

#include "WallDist.h"
#include <memory>
#include "CmxTask.h"
#include "TaskState.h"
#include "NsCtrl.h"
#include "SimuDef.h"
#include "SolverDef.h"
#include "SolverState.h"
#include "Parallel.h"
#include "Zone.h"
#include "ZoneState.h"
#include "UnsGrid.h"
#include "FaceTopo.h"
#include "CellMesh.h"
#include "FaceMesh.h"
#include "NodeMesh.h"
#include "BcRecord.h"
#include "Boundary.h"
#include "ActionState.h"
#include "DataBook.h"
#include "DataBaseIO.h"
#include "InterFace.h"
#include "GteVector.h"
#include "GteSegment.h"
#include "GteDCPQuery.h"
#include "GteDistPointSegment.h"
#include "GteTriangle.h"
#include "GteDistPointTriangleExact.h"
#include "HXMath.h"
#include "LogFile.h"
#include <utility>

BeginNameSpace( ONEFLOW )

namespace
{
    // Process-wide storage for the FILL -> CALC wall-distance pipeline.
    // unique_ptr guarantees a single owner and exception-safe cleanup.
    std::unique_ptr< WallStructure > g_wallStructure;
}

WallStructure * GetWallStructure() noexcept
{
    return g_wallStructure.get();
}

void ResetWallStructure()
{
    g_wallStructure.reset();
}

void EnsureWallStructure()
{
    if ( ! g_wallStructure )
    {
        g_wallStructure = std::make_unique< WallStructure >();
    }
}

void FreeWallStruct()
{
    ResetWallStructure();
}

void SetWallTask()
{
    REGISTER_DATA_CLASS( FillWallStructTask );
    REGISTER_DATA_CLASS( FillWallStruct );
    REGISTER_DATA_CLASS( CalcWallDist );
}

void FillWallStructTask( StringField & /*data*/ )
{
    // TaskState still expects a raw Task*; ownership is transferred there.
    TaskState::createdTask = std::make_unique<CFillWallStructTaskImp>();
}

void FillWallStruct( StringField & /*data*/ )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    const int nBFaces = grid->GetFaceTopo().bcManager->bcRecord->GetNBFace();
    BcRecord * bcRecord = grid->GetFaceTopo().bcManager->bcRecord.get();

    const int nWallFace = bcRecord->CalcNumWallFace();

    RealField & xfc = grid->GetFaceMesh().xfc;
    RealField & yfc = grid->GetFaceMesh().yfc;
    RealField & zfc = grid->GetFaceMesh().zfc;

    RealField & x = grid->nodeMesh->xN;
    RealField & y = grid->nodeMesh->yN;
    RealField & z = grid->nodeMesh->zN;

    ActionState::dataBook->MoveToBegin();
    HXWrite( ActionState::dataBook, nWallFace );

    if ( nWallFace <= 0 )
    {
        return;
    }

    WallStructure::PointField fc;
    WallStructure::PointLink  fv;

    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        const int bcType = bcRecord->bcType[ iFace ];
        const int nNodes = static_cast< int >( grid->GetFaceTopo().faces[ iFace ].size() );

        if ( bcType != BC::SOLID_SURFACE )
        {
            continue;
        }

        WallStructure::PointField simpleFace;
        simpleFace.reserve( static_cast< size_t >( nNodes ) );

        for ( int iNode = 0; iNode < nNodes; ++ iNode )
        {
            const int iPoint = grid->GetFaceTopo().faces[ iFace ][ iNode ];
            simpleFace.emplace_back( x[ iPoint ], y[ iPoint ], z[ iPoint ] );
        }
        fv.push_back( std::move( simpleFace ) );
        fc.emplace_back( xfc[ iFace ], yfc[ iFace ], zfc[ iFace ] );
    }

    HXWrite( ActionState::dataBook, fv );
    HXWrite( ActionState::dataBook, fc );
}

void CalcWallDist( StringField & /*data*/ )
{
    WallStructure * ws = GetWallStructure();
    if ( ! ws )
    {
        return;
    }

    UnsGrid * grid = Zone::GetUnsGrid();
    RealField & dist = grid->GetCellMesh().dist;
    const int nCells = grid->nCells;

    dist = LARGE;

    RealField & xcc = grid->GetCellMesh().xcc;
    RealField & ycc = grid->GetCellMesh().ycc;
    RealField & zcc = grid->GetCellMesh().zcc;

    std::cout << "zone " << grid->id << std::endl;

    WallStructure::PointField & fc = ws->fc;
    WallStructure::PointLink  & fv = ws->fv;
    const int nWFace = static_cast< int >( fc.size() );

    for ( int cId = 0; cId < nCells; ++ cId )
    {
        if ( cId % 10000 == 0 )
        {
            std::cout << " pid = " << Parallel::pid << " Zone = " << grid->id
                      << " cid = " << cId << " nCells = " << nCells
                      << " nWFace = " << nWFace << std::endl;
        }

        WallStructure::PointType ccp( xcc[ cId ], ycc[ cId ], zcc[ cId ] );

        for ( int iWFace = 0; iWFace < nWFace; ++ iWFace )
        {
            WallStructure::PointField & fvList = fv[ iWFace ];
            const Real wdst = CalcPoint2FaceDist( ccp, fvList );
            if ( dist[ cId ] > wdst )
            {
                dist[ cId ] = wdst;
            }
        }
    }

    for ( int cId = 0; cId < nCells; ++ cId )
    {
        dist[ cId ] = sqrt( dist[ cId ] );
    }
}

void CFillWallStructTaskImp::Run()
{
    ActionState::dataBook = this->dataBook.get();
    this->Create();

    for ( int zId = 0; zId < ZoneState::nZones; ++ zId )
    {
        ZoneState::zid = zId;

        if ( Parallel::pid == ZoneState::pid[ zId ] )
        {
            this->action();
        }

        HXBcast( ActionState::dataBook, ZoneState::pid[ zId ] );
        FillWall();
    }
}

void CFillWallStructTaskImp::Create()
{
    // Fresh storage for this aggregation pass.
    g_wallStructure = std::make_unique< WallStructure >();
}

void CFillWallStructTaskImp::FillWall()
{
    WallStructure * ws = GetWallStructure();
    if ( ! ws )
    {
        return;
    }

    ActionState::dataBook->MoveToBegin();
    int nSolidCells = 0;
    HXRead( ActionState::dataBook, nSolidCells );

    WallStructure::PointField fcTmp;
    WallStructure::PointLink  fvTmp;

    fvTmp.resize( static_cast< size_t >( nSolidCells ) );
    HXRead( ActionState::dataBook, fvTmp );

    fcTmp.resize( static_cast< size_t >( nSolidCells ) );
    HXRead( ActionState::dataBook, fcTmp );

    for ( int cId = 0; cId < nSolidCells; ++ cId )
    {
        ws->fc.push_back( fcTmp[ cId ] );
        ws->fv.push_back( std::move( fvTmp[ cId ] ) );
    }
}

Real CalcPoint2FaceDist( WallStructure::PointType node,
                         WallStructure::PointField & fvList )
{
    using namespace gte;

    Vector< 3, Real > point0;
    point0[ 0 ] = node.x;
    point0[ 1 ] = node.y;
    point0[ 2 ] = node.z;

    Vector< 3, Real > point1;
    Vector< 3, Real > point2;
    Vector< 3, Real > point3;

    const int nVertex = static_cast< int >( fvList.size() );

    if ( nVertex <= 2 )
    {
        point1[ 0 ] = fvList[ 0 ].x;
        point1[ 1 ] = fvList[ 0 ].y;
        point1[ 2 ] = fvList[ 0 ].z;

        point2[ 0 ] = fvList[ 1 ].x;
        point2[ 1 ] = fvList[ 1 ].y;
        point2[ 2 ] = fvList[ 1 ].z;

        Segment< 3, Real > segment( point1, point2 );
        using SuperLine = DCPQuery< Real, Vector< 3, Real >, Segment< 3, Real > >;
        SuperLine query;
        return query( point0, segment ).sqrDistance;
    }

    if ( nVertex == 3 )
    {
        point1[ 0 ] = fvList[ 0 ].x;
        point1[ 1 ] = fvList[ 0 ].y;
        point1[ 2 ] = fvList[ 0 ].z;

        point2[ 0 ] = fvList[ 1 ].x;
        point2[ 1 ] = fvList[ 1 ].y;
        point2[ 2 ] = fvList[ 1 ].z;

        point3[ 0 ] = fvList[ 2 ].x;
        point3[ 1 ] = fvList[ 2 ].y;
        point3[ 2 ] = fvList[ 2 ].z;

        Triangle< 3, Real > triangle( point1, point2, point3 );
        DistancePointTriangleExact< 3, Real > query;
        return query( point0, triangle ).sqrDistance;
    }

    Real xCenter = 0;
    Real yCenter = 0;
    Real zCenter = 0;
    for ( int iv = 0; iv < nVertex; ++ iv )
    {
        xCenter += fvList[ iv ].x;
        yCenter += fvList[ iv ].y;
        zCenter += fvList[ iv ].z;
    }
    const Real coef = 1.0 / nVertex;
    xCenter *= coef;
    yCenter *= coef;
    zCenter *= coef;

    point3[ 0 ] = xCenter;
    point3[ 1 ] = yCenter;
    point3[ 2 ] = zCenter;

    Real dist = LARGE;
    for ( int iv = 0; iv < nVertex; ++ iv )
    {
        const int p0 = iv;
        const int p1 = ( iv + 1 ) % nVertex;

        point1[ 0 ] = fvList[ p0 ].x;
        point1[ 1 ] = fvList[ p0 ].y;
        point1[ 2 ] = fvList[ p0 ].z;

        point2[ 0 ] = fvList[ p1 ].x;
        point2[ 1 ] = fvList[ p1 ].y;
        point2[ 2 ] = fvList[ p1 ].z;

        Triangle< 3, Real > triangle( point1, point2, point3 );
        DistancePointTriangleExact< 3, Real > query;
        dist = MIN( dist, query( point0, triangle ).sqrDistance );
    }
    return dist;
}

EndNameSpace

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

#include "WallDist.h"
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

BeginNameSpace( ONEFLOW )

WallStructure * wallstruct = 0;

void FreeWallStruct()
{
    delete wallstruct;
}

void SetWallTask()
{
    REGISTER_DATA_CLASS( FillWallStructTask  );
    REGISTER_DATA_CLASS( FillWallStruct  );
    REGISTER_DATA_CLASS( CalcWallDist );
}

void FillWallStructTask( StringField & data )
{
    CFillWallStructTaskImp * task = new CFillWallStructTaskImp();
    TaskState::task = task;
}

void FillWallStruct( StringField & data )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    int nBFace = grid->faceTopo->bcManager->bcRecord->GetNBFace();
    BcRecord * bcRecord = grid->faceTopo->bcManager->bcRecord;

    int nWallFace = bcRecord->CalcNumWallFace();

    RealField & xfc = grid->faceMesh->xfc;
    RealField & yfc = grid->faceMesh->yfc;
    RealField & zfc = grid->faceMesh->zfc;

    RealField & x = grid->nodeMesh->xN;
    RealField & y = grid->nodeMesh->yN;
    RealField & z = grid->nodeMesh->zN;

    ActionState::dataBook->MoveToBegin();

    HXWrite( ActionState::dataBook, nWallFace );

    if ( nWallFace  <= 0 )
    {
        return;
    }

    WallStructure::PointField fc;
    WallStructure::PointLink  fv;

    for ( int iFace = 0; iFace < nBFace; ++ iFace )
    {
        int bcType = bcRecord->bcType[ iFace ];
        int bcRegion = bcRecord->bcRegion[ iFace ];
        int nNode = grid->faceTopo->f2n[ iFace ].size();

        if ( bcType == BC::SOLID_SURFACE )
        {
            WallStructure::PointField simpleFace;

            for ( int iNode = 0; iNode < nNode; ++ iNode )
            {
                int iPoint = grid->faceTopo->f2n[ iFace ][ iNode ];
                Real x0 = x[ iPoint ];
                Real y0 = y[ iPoint ];
                Real z0 = z[ iPoint ];
                WallStructure::PointType simplePoint( x0, y0, z0 );
                simpleFace.push_back( simplePoint );
            }
            fv.push_back( simpleFace );

            Real xxfc = xfc[ iFace ];
            Real yyfc = yfc[ iFace ];
            Real zzfc = zfc[ iFace ];
            WallStructure::PointType centerPoint ( xxfc, yyfc, zzfc );
            fc.push_back( centerPoint );
        }
    }

    HXWrite( ActionState::dataBook, fv );
    HXWrite( ActionState::dataBook, fc );
}

void CalcWallDist( StringField & data )
{
    UnsGrid * grid = Zone::GetUnsGrid();
    int nBFace = grid->faceTopo->bcManager->bcRecord->GetNBFace();
    BcRecord * bcRecord = grid->faceTopo->bcManager->bcRecord;

    int nWallFace = bcRecord->CalcNumWallFace();

    RealField & dist = grid->cellMesh->dist;

    int nCell = grid->nCell;

    dist = LARGE;

    RealField & xcc = grid->cellMesh->xcc;
    RealField & ycc = grid->cellMesh->ycc;
    RealField & zcc = grid->cellMesh->zcc;

    cout << "zone " << grid->id << endl;

    WallStructure::PointField & fc = wallstruct->fc;
    WallStructure::PointLink  & fv = wallstruct->fv;

    int nWFace = wallstruct->fc.size();

    for ( int cId = 0; cId < nCell; ++ cId )
    {
        if ( cId % 10000 == 0 )
        {
            cout << " pid = " << Parallel::pid << " Zone = " << grid->id;
            cout << " cid = " << cId << " " << "nCell = " << nCell;
            cout << " nWFace = " << nWFace << endl;
        }

        Real xc = xcc[ cId ];
        Real yc = ycc[ cId ];
        Real zc = zcc[ cId ];

        WallStructure::PointType ccp( xc, yc, zc );

        for ( int iWFace = 0; iWFace < nWFace; ++ iWFace )
        {
            WallStructure::PointField  & fvList = fv[ iWFace ];

            Real wdst = CalcPoint2FaceDist( ccp, fvList );

            if ( wdst < 1.0e-30 )
            {
            }

            if ( dist[ cId ] > wdst )
            {
                dist[ cId ] = wdst;
            }
        }        
    }


    for ( int cId = 0; cId < nCell; ++ cId )
    {
        dist[ cId ] = sqrt( dist[ cId ] );
    }
}


CFillWallStructTaskImp::CFillWallStructTaskImp()
{
    ;
}

CFillWallStructTaskImp::~CFillWallStructTaskImp()
{
    ;
}

void CFillWallStructTaskImp::Run()
{
    ActionState::dataBook = this->dataBook;
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
    wallstruct = new WallStructure();
}

void CFillWallStructTaskImp::FillWall()
{
    ActionState::dataBook->MoveToBegin();
    int nSolidCell;
    HXRead( ActionState::dataBook, nSolidCell );

    WallStructure::PointField fcTmp;
    WallStructure::PointLink  fvTmp;

    fvTmp.resize( nSolidCell );
    HXRead( ActionState::dataBook, fvTmp );

    fcTmp.resize( nSolidCell );
    HXRead( ActionState::dataBook, fcTmp );

    for ( int cId = 0; cId < nSolidCell; ++ cId )
    {
        wallstruct->fc.push_back( fcTmp[ cId ] );
        wallstruct->fv.push_back( fvTmp[ cId ] );
    }
}

Real CalcPoint2FaceDist( WallStructure::PointType node, WallStructure::PointField & fvList )
{
    using namespace gte;
    Vector< 3, Real > point0;

    point0[ 0 ] = node.x;
    point0[ 1 ] = node.y;
    point0[ 2 ] = node.z;

    Vector< 3, Real > point1;
    Vector< 3, Real > point2;
    Vector< 3, Real > point3;

    int nVertex = fvList.size();

    if ( nVertex <= 2 )
    {
        point1[ 0 ] = fvList[ 0 ].x;
        point1[ 1 ] = fvList[ 0 ].y;
        point1[ 2 ] = fvList[ 0 ].z;

        point2[ 0 ] = fvList[ 1 ].x;
        point2[ 1 ] = fvList[ 1 ].y;
        point2[ 2 ] = fvList[ 1 ].z;

        Segment< 3, Real > segment( point1, point2 );

        typedef DCPQuery<Real, Vector< 3, Real >, Segment< 3, Real > > SuperLine;

        SuperLine b;

        SuperLine::Result result = b( point0, segment );
        return result.sqrDistance;
    }
    else if ( nVertex == 3 )
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

        DistancePointTriangleExact< 3, Real > a;
        DistancePointTriangleExact< 3, Real >::Result r = a( point0, triangle );
        return r.sqrDistance;
    }
    else
    {
        Real xCenter = 0;
        Real yCenter = 0;
        Real zCenter = 0;
        for ( int iv = 0; iv < nVertex; ++ iv )
        {
            xCenter += fvList[ iv ].x;
            yCenter += fvList[ iv ].y;
            zCenter += fvList[ iv ].z;
        }
        Real coef = 1.0 / nVertex;
        xCenter *= coef;
        yCenter *= coef;
        zCenter *= coef;

        point3[ 0 ] = xCenter;
        point3[ 1 ] = yCenter;
        point3[ 2 ] = zCenter;

        Real dist = LARGE;

        for ( int iv = 0; iv < nVertex; ++ iv )
        {
            int p0 = iv;
            int p1 = ( iv + 1 ) % nVertex;

            point1[ 0 ] = fvList[ p0 ].x;
            point1[ 1 ] = fvList[ p0 ].y;
            point1[ 2 ] = fvList[ p0 ].z;

            point2[ 0 ] = fvList[ p1 ].x;
            point2[ 1 ] = fvList[ p1 ].y;
            point2[ 2 ] = fvList[ p1 ].z;

            Triangle< 3, Real > triangle( point1, point2, point3 );

            DistancePointTriangleExact< 3, Real > a;
            DistancePointTriangleExact< 3, Real >::Result r = a( point0, triangle );
            dist = MIN( dist, r.sqrDistance );
        }
        return dist;
    }
}

EndNameSpace
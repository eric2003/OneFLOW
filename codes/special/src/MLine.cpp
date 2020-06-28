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

#include "MLine.h"
#include "CurveMesh.h"
#include "BlockFaceSolver.h"
#include "SimpleDomain.h"
#include "LineMachine.h"
#include "DomainMachine.h"
#include "BlkMesh.h"
#include "SDomain.h"
#include "HXPointer.h"
#include "CurveInfo.h"
#include "SegmentCtrl.h"
#include "BlockElem.h"
#include "HXCgns.h"
#include "Dimension.h"
#include <algorithm>
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )


SLine::SLine()
{
    ;
}

SLine::~SLine()
{
}

void SLine::Alloc()
{
    this->x1d.resize( ni );
    this->y1d.resize( ni );
    this->z1d.resize( ni );
}

void SLine::CopyMesh()
{
    CurveMesh * curveMesh = line_Machine.GetCurveMesh( line_id );
    int p1 = curveMesh->curveInfo->p1;
    int p2 = curveMesh->curveInfo->p2;
    this->ctrlpoints.push_back( p1 );
    this->ctrlpoints.push_back( p2 );
    for ( int i = 1; i <= ni; ++ i )
    {
        int i0 = i - 1;
        PointType * pt = curveMesh->ptList[ i0 ];
        this->x1d[ i0 ] = pt->x;
        this->y1d[ i0 ] = pt->y;
        this->z1d[ i0 ] = pt->z;
    }
    ;
}

void SLine::ConstructCtrlPoints()
{
    IntField & pointIdList = blkFaceSolver.myFaceSolver.GetLine( line_id );
    this->ctrlpoints = pointIdList;
}

void SLine::SetDomainBcMesh( SDomain * sDomain )
{
    RealField2D & x2d = sDomain->x2d;
    RealField2D & y2d = sDomain->y2d;
    RealField2D & z2d = sDomain->z2d;

    int line_id = this->line_id - 1;
    SLine * sLine = blkFaceSolver.myFaceSolver.slineList[ line_id ];
    ni = sLine->ni;
    RealField & x1d = sLine->x1d;
    RealField & y1d = sLine->y1d;
    RealField & z1d = sLine->z1d;

    int p1 = this->ctrlpoints[ 0 ];
    int p2 = this->ctrlpoints[ 1 ];

    int q1 = sLine->ctrlpoints[ 0 ];
    int q2 = sLine->ctrlpoints[ 1 ];

    int flag = 1;
    if ( p1 != q1 ) flag = -1;

    CoorMap::iterator it1 = sDomain->localCoorMap->find( p1 );
    CoorMap::iterator it2 = sDomain->localCoorMap->find( p2 );

    CalcCoor & c1 = it1->second;
    CalcCoor & c2 = it2->second;

    CalcCoor d1;
    d1.i = c2.i - c1.i;
    d1.j = c2.j - c1.j;
    d1.k = c2.k - c1.k;

    GetUnitDir( d1 );

    for ( int i = 1; i <= ni; ++ i )
    {
        CalcCoor ci;
        int i0 = i - 1;
        ci.i = c1.i + d1.i * i0;
        ci.j = c1.j + d1.j * i0;
        ci.k = c1.k + d1.k * i0;

        int ii = ci.i - 1;
        int jj = ci.j - 1;
        int kk = ci.k - 1;

        if ( flag == -1 )
        {
            i0 = ni - 1 - i0;
        }
        Real xm = x1d[ i0 ];
        Real ym = y1d[ i0 ];
        Real zm = z1d[ i0 ];
        x2d[ ii ][ jj ] = xm;
        y2d[ ii ][ jj ] = ym;
        z2d[ ii ][ jj ] = zm;
    }

    int kkk = 1;
}

MLine::MLine()
{
    ;
}

MLine::~MLine()
{
    DeletePointer( slineList );
}

void MLine::ConstructSLineCtrlPoint()
{
    int nSline = this->slineList.size();
    for ( int iSLine = 0; iSLine < nSline; ++ iSLine )
    {
        SLine * sLine = this->slineList[ iSLine ];
        sLine->ConstructCtrlPoints();
    }
}

void MLine::ConstructCtrlPoint()
{
    map< int, IntSet >::iterator iter;
    for ( iter = this->pointToDomainMap.begin(); iter != this->pointToDomainMap.end(); ++ iter )
    {
        if ( iter->second.size() == 1 )
        {
            this->candidate_ctrlpoints.push_back( iter->first );
        }
    }
}

void MLine::ConstructPointToPointMap()
{
    this->ConstructPointToPointMap( this->pointToPointMap );
}

void MLine::ConstructPointToPointMap( map< int, IntSet > & pointToPointMap )
{
    MLine * mLine = this;
    LinkField pointIdLink;
    GetPointIdLink( mLine->lineList, pointIdLink );

    ONEFLOW::ConstructPointToPointMap( pointIdLink, pointToPointMap );
}

void MLine::ConstructPointToDomainMap()
{
    for ( int iLine = 0; iLine < lineList.size(); ++ iLine )
    {
        int line_id = lineList[ iLine ];
        IntField & pointIdList = blkFaceSolver.myFaceSolver.GetLine( line_id );

        ConstructIntList2Map( line_id, pointIdList, pointToDomainMap );
    }
}

void MLine::ConstructPointToDomainMap( int domain_id, map< int, IntSet > & pointToDomainMap )
{
    MLine * mLine = this;
    LinkField pointIdLink;
    GetPointIdLink( mLine->lineList, pointIdLink );

    ONEFLOW::ConstructPointToDomainMap( domain_id, pointIdLink, pointToDomainMap );
}

void MLine::ConstructLineToDomainMap()
{
    for ( int iLine = 0; iLine < lineList.size(); ++ iLine )
    {
        int line_id = lineList[ iLine ];

        ConstructInt2Map( line_id, line_id, this->lineToDomainMap );
    }
}

void MLine::ConstructLineToDomainMap( int domain_id, map< int, IntSet > & lineToDomainMap )
{
    MLine * mLine = this;
    ONEFLOW::ConstructLineToDomainMap( domain_id, mLine->lineList, lineToDomainMap );
}

void MLine::CalcCoor( CoorMap * localCoorMap )
{
    int openLine = 0;
    this->CalcBcCoor( localCoorMap, openLine );
}

void MLine::ConstructDomainTopo()
{
    this->ConstructLineToDomainMap();
    this->ConstructPointToDomainMap();
    this->ConstructPointToPointMap();
    this->ConstructBcPoint();
    this->ConstructCtrlPoint();
}

void MLine::AddSubLine( int line_id )
{
    this->lineList.push_back( line_id );
    SLine * sLine = new SLine();
    sLine->line_id = line_id;
    this->slineList.push_back( sLine );
}

void MLine::SetDomainBcMesh( SDomain * sDomain )
{
    int nSLine = slineList.size();
    for ( int iSLine = 0; iSLine < nSLine; ++ iSLine )
    {
        SLine * sLine = this->slineList[ iSLine ];
        sLine->SetDomainBcMesh( sDomain );
    }
}

void MLine::CreateInpFaceList( HXVector< Face2D * > &facelist )
{
    for ( int iSLine = 0; iSLine < this->slineList.size(); ++ iSLine )
    {
        SLine * sLine = this->slineList[ iSLine ];
        Face2D * face2d = new Face2D();
        face2d->face_id = sLine->line_id;
        face2d->Set1DRegion( sLine->ctrlpoints );
        BlkF2C & face_struct = blkFaceSolver.myFaceSolver.face2Block[ face2d->face_id ];
        face2d->bcType = face_struct.bctype;
        face2d->CalcStEd( coorMap );
        facelist.push_back( face2d );
    }
}



EndNameSpace
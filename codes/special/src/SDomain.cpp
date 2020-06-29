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

#include "BlockFaceSolver.h"
#include "BlkMesh.h"
#include "MLine.h"
#include "Transfinite.h"
#include "SimpleDomain.h"
#include "LineMachine.h"
#include "DomainMachine.h"
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

SDomain::SDomain()
{
    int nMLine = 4;
    for ( int iMLine = 0; iMLine < nMLine; ++ iMLine )
    {
        MLine * mLine = new MLine();
        mLine->pos = iMLine;
        mLineList.push_back( mLine );
    }
    localCoorMap = new CoorMap();
}

SDomain::~SDomain()
{
    DeletePointer( mLineList );
    delete localCoorMap;
}

void SDomain::Alloc()
{
    AllocateVector( x2d, ni, nj );
    AllocateVector( y2d, ni, nj );
    AllocateVector( z2d, ni, nj );
}

void SDomain::SetDomain( int fid, IntField & lineList, IntField & posList )
{
    this->domain_id = fid;

    for ( int iLine = 0; iLine < lineList.size(); ++ iLine )
    {
        int line_id = lineList[ iLine ];
        int pos = posList[ iLine ];
        MLine * mLine = mLineList[ pos ];
        mLine->AddSubLine( line_id );
    }
}

void SDomain::ConstructSDomainCtrlPoint()
{
    int nMLine = mLineList.size();
    for ( int iMLine = 0; iMLine < nMLine; ++ iMLine )
    {
        MLine * mLine = mLineList[ iMLine ];
        mLine->ConstructPointToDomainMap();
        mLine->ConstructCtrlPoint();
        mLine->ConstructSLineCtrlPoint();
    }

    for ( int iMLine = 0; iMLine < nMLine; ++ iMLine )
    {
        int iMLine1 = iMLine - 1;
        if ( iMLine1 < 0 ) iMLine1 = nMLine - 1;
        MLine * mLine1 = mLineList[ iMLine1 ];
        MLine * mLine2 = mLineList[ iMLine  ];
        int pt;
        GetCommonPoint( mLine1, mLine2, pt );
        this->ctrlpoints.push_back( pt );
    }

    for ( int iMLine = 0; iMLine < nMLine; ++ iMLine )
    {
        MLine * mLine = mLineList[ iMLine  ];
        int p1 = this->ctrlpoints[ iMLine ];
        int p2 = this->ctrlpoints[ ( iMLine + 1 ) % nMLine ];
        mLine->ctrlpoints.push_back( p1 );
        mLine->ctrlpoints.push_back( p2 );
    }
}

void SDomain::GetCommonPoint( MLine * mLine1, MLine * mLine2, int & pt )
{
    int p1 = mLine1->candidate_ctrlpoints[ 0 ];
    int p2 = mLine1->candidate_ctrlpoints[ 1 ];

    int q1 = mLine2->candidate_ctrlpoints[ 0 ];
    int q2 = mLine2->candidate_ctrlpoints[ 1 ];

    if ( p1 == q1 )
    {
        pt = p1;
    }
    else if ( p1 == q2 )
    {
        pt = p1;
    }
    else
    {
        pt = p2;
    }
}

void SDomain::SetRemainingCtrlPoint( IntField & idxList )
{
    int nCtrlPoint = this->ctrlpoints.size();
    int ip = -1;
    for ( int i = 0; i < nCtrlPoint; ++ i )
    {
        if ( idxList[ i ] != 1 )
        {
            ip = i;
            break;
        }
    }
    int ip1 = ( ip + 1 ) % nCtrlPoint;
    int ip2 = ( ip + 2 ) % nCtrlPoint;
    int ip3 = ( ip + 3 ) % nCtrlPoint;

    CoorMap::iterator it1 = this->coorMap->find( ip1 );
    CoorMap::iterator it2 = this->coorMap->find( ip2 );
    CoorMap::iterator it3 = this->coorMap->find( ip3 );

    int mi = it1->second.i + it3->second.i - it2->second.i;
    int mj = it1->second.j + it3->second.j - it2->second.j;
    int mk = it1->second.k + it3->second.k - it2->second.k;

    int pt = this->ctrlpoints[ ip ];

    CalcCoor coor;
    coor.SetCoor( mi, mj, mk );
    coorMap->insert( pair<int, CalcCoor>( pt, coor ) );
}


bool SDomain::CalcSingleDomainCoor()
{
    int nCtrlPoint = this->ctrlpoints.size();
    int nCount = 0;
    IntField idx( nCtrlPoint, 0 );
    for ( int i = 0; i < nCtrlPoint; ++ i )
    {
        int ip = this->ctrlpoints[ i ];

        CoorMap::iterator iter = this->coorMap->find( ip );
        if ( iter != this->coorMap->end() )
        {
            idx[ i ] = 1;
            ++ nCount;
        }
    }

    if ( nCount < 3 ) return false;

    if ( nCount == 3 )
    {
        this->SetRemainingCtrlPoint( idx );
        return false;
    }
    else if ( nCount == 4 )
    {
        this->CalcDim2D();
        int closedLine = 1;
        this->CalcBcCoor( this->coorMap, closedLine );
        return true;
    }
    return true;
}

void SDomain::ConstructPointToPointMap()
{
    this->ConstructPointToPointMap( this->pointToPointMap );
}

void SDomain::ConstructPointToPointMap( map< int, IntSet > & pointToPointMap )
{
    for ( int iMLine = 0; iMLine < mLineList.size(); ++ iMLine )
    {
        MLine * mLine = mLineList[ iMLine ];
        LinkField pointIdLink;
        this->GetPointIdLink( mLine->lineList, pointIdLink );

        ONEFLOW::ConstructPointToPointMap( pointIdLink, pointToPointMap );
    }
}

void SDomain::ConstructPointToDomainMap()
{
    this->ConstructPointToDomainMap( this->pointToDomainMap );
}

void SDomain::GetPointIdLink( IntField & lineList, LinkField & pointIdLink )
{
    for ( int iLine = 0; iLine < lineList.size(); ++ iLine )
    {
        int line_id = lineList[ iLine ] - 1;
        IntField & pointIdList = blkFaceSolver.myFaceSolver.lineList[ line_id ];
        pointIdLink.push_back( pointIdList );
    }
}

void SDomain::ConstructPointToDomainMap( map< int, IntSet > & pointToDomainMap )
{
    for ( int iMLine = 0; iMLine < mLineList.size(); ++ iMLine )
    {
        MLine * mLine = mLineList[ iMLine ];
        LinkField pointIdLink;
        this->GetPointIdLink( mLine->lineList, pointIdLink );

        ONEFLOW::ConstructPointToDomainMap( this->domain_id, pointIdLink, pointToDomainMap );
    }
}

void SDomain::ConstructLineToDomainMap()
{
    this->ConstructLineToDomainMap( this->lineToDomainMap );
}

void SDomain::ConstructLineToDomainMap( map< int, IntSet > & lineToDomainMap )
{
    for ( int iMLine = 0; iMLine < mLineList.size(); ++ iMLine )
    {
        MLine * mLine = mLineList[ iMLine ];
        ONEFLOW::ConstructLineToDomainMap( this->domain_id, mLine->lineList, lineToDomainMap );
    }
}

void SDomain::ConstructDomainTopo()
{
    this->ConstructLineToDomainMap();
    this->ConstructPointToDomainMap();
    this->ConstructPointToPointMap();
    this->ConstructBcPoint();
    this->ConstructCtrlPoint();
}

void SDomain::Add( IntField &iList, IntField &jList, IntField &kList, int i, int j, int k )
{
    iList.push_back( i );
    jList.push_back( j );
    kList.push_back( k );
}

void SDomain::ConstructLocalTopoAsBlk2D()
{
    IntField iList, jList, kList;
    Add( iList, jList, kList, 1, 1, 1 );
    Add( iList, jList, kList, ni, 1, 1 );
    Add( iList, jList, kList, ni, nj, 1 );
    Add( iList, jList, kList, 1, nj, 1 );

    for ( int iPoint = 0; iPoint < iList.size(); ++ iPoint )
    {
        int pt = this->ctrlpoints[ iPoint ];
        int i = iList[ iPoint ];
        int j = jList[ iPoint ];
        int k = kList[ iPoint ];
        CalcCoor c;
        c.SetCoor( i, j, k );
        this->localCoorMap->insert( pair<int, CalcCoor>( pt, c ) );
    }

    int nMLine = mLineList.size();
    for ( int iMLine = 0; iMLine < nMLine; ++ iMLine )
    {
        MLine * mLine = mLineList[ iMLine ];
        mLine->ConstructDomainTopo();
        mLine->CalcDim1D();
        mLine->CalcCoor( this->localCoorMap );
    }
}

void SDomain::SetBlkBcMesh( Block3D * blk3d )
{
    RealField3D & x3d = blk3d->x3d;
    RealField3D & y3d = blk3d->y3d;
    RealField3D & z3d = blk3d->z3d;

    SDomain * sDomain = blkFaceSolver.myFaceSolver.sDomainList[ this->domain_id ];

    RealField2D & x2d = sDomain->x2d;
    RealField2D & y2d = sDomain->y2d;
    RealField2D & z2d = sDomain->z2d;

    int p1 = this->ctrlpoints[ 0 ];
    int p2 = this->ctrlpoints[ 1 ];
    int p3 = this->ctrlpoints[ 2 ];
    int p4 = this->ctrlpoints[ 3 ];

    int q1 = sDomain->ctrlpoints[ 0 ];
    int q2 = sDomain->ctrlpoints[ 1 ];
    int q3 = sDomain->ctrlpoints[ 2 ];
    int q4 = sDomain->ctrlpoints[ 3 ];

    CoorMap::iterator it1 = coorMap->find( p1 );
    CoorMap::iterator it2 = coorMap->find( p2 );
    CoorMap::iterator it3 = coorMap->find( p3 );
    CoorMap::iterator it4 = coorMap->find( p4 );

    CalcCoor & c1 = it1->second;
    CalcCoor & c2 = it2->second;
    CalcCoor & c3 = it3->second;
    CalcCoor & c4 = it4->second;

    CalcCoor d1;
    d1.i = c2.i - c1.i;
    d1.j = c2.j - c1.j;
    d1.k = c2.k - c1.k;

    CalcCoor d2;
    d2.i = c4.i - c1.i;
    d2.j = c4.j - c1.j;
    d2.k = c4.k - c1.k;

    GetUnitDir( d1 );
    GetUnitDir( d2 );

    int ii, jj, kk;
    int i0, j0;

    for ( int j = 1; j <= nj; ++ j )
    {
        CalcCoor cj;
        j0 = j - 1;
        cj.i = c1.i + d2.i * j0;
        cj.j = c1.j + d2.j * j0;
        cj.k = c1.k + d2.k * j0;

        for ( int i = 1; i <= ni; ++ i )
        {
            i0 = i - 1;
            CalcCoor ci;
            ci.i = c1.i + d1.i * i0;
            ci.j = c1.j + d1.j * i0;
            ci.k = c1.k + d1.k * i0;
            ii = ci.i + cj.i - c1.i - 1;
            jj = ci.j + cj.j - c1.j - 1;
            kk = ci.k + cj.k - c1.k - 1;
            //cout << " i,j = " << i << " " << j << " ni, nj = " << ni << " " << nj << "\n";
            //cout << " ii,jj,kk = " << ii << " " << jj << " " << kk << "\n";
            Real xm = x2d[ i0 ][ j0 ];
            Real ym = y2d[ i0 ][ j0 ];
            Real zm = z2d[ i0 ][ j0 ];
            x3d[ ii ][ jj ][ kk ] = xm;
            y3d[ ii ][ jj ][ kk ] = ym;
            z3d[ ii ][ jj ][ kk ] = zm;
        }
    }
    int kkk = 1;
}

void SDomain::SetDomainBcMesh()
{
    int nMLine = mLineList.size();
    for ( int iMLine = 0; iMLine < nMLine; ++ iMLine )
    {
        MLine * mLine = mLineList[ iMLine ];
        mLine->SetDomainBcMesh( this );
    }
}

void SDomain::GenerateSDomainMesh()
{
    TransfiniteInterpolation( x2d, ni, nj );
    TransfiniteInterpolation( y2d, ni, nj );
    TransfiniteInterpolation( z2d, ni, nj );
}

void SDomain::GenerateSDomainMesh( fstream & file )
{
    TransfiniteInterpolation( x2d, ni, nj );
    TransfiniteInterpolation( y2d, ni, nj );
    TransfiniteInterpolation( z2d, ni, nj );

    file << " VARIABLES = \"X\", \"Y\", \"Z\" \n";
    file << " ZONE I = " << ni << ", J = " << nj << " F = POINT \n";
    for ( int j = 0; j < nj; ++ j )
    {
        for ( int i = 0; i < ni; ++ i )
        {
            file << x2d[ i ][ j ] << " " << y2d[ i ][ j ] << " " << z2d[ i ][ j ] << "\n";
        }
    }
}



EndNameSpace
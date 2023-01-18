/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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
#include "SimpleDomain.h"
#include "BlockFaceSolver.h"
#include "HXStd.h"
#include "LineMachine.h"
#include "DomainMachine.h"
#include "CurveInfo.h"
#include "SegmentCtrl.h"
#include "BlockElem.h"
#include "HXCgns.h"
#include "Dimension.h"
#include <algorithm>
#include <iostream>


BeginNameSpace( ONEFLOW )

BlkF2C::BlkF2C()
{
    ;
}

BlkF2C::~BlkF2C()
{
    ;
}

Face2D::Face2D()
{
    this->t = 0;
}

Face2D::~Face2D()
{
    delete this->t;
}

void Face2D::CalcRegion()
{
    if ( this->bcType != -1 )
    {
        p1 = st;
        p2 = ed;
        this->st.i = MIN( p1.i, p2.i );
        this->st.j = MIN( p1.j, p2.j );
        this->st.k = MIN( p1.k, p2.k );

        this->ed.i = MAX( p1.i, p2.i );
        this->ed.j = MAX( p1.j, p2.j );
        this->ed.k = MAX( p1.k, p2.k );
    }
}

void Face2D::CalcStEd( CoorMap * coorMap )
{
    CoorMap::iterator c1, c2;
    int p1 = this->ctrlpoints[ 0 ];
    int p3 = this->ctrlpoints[ 2 ];
    c1 = coorMap->find( p1 );
    c2 = coorMap->find( p3 );
    st = c1->second;
    ed = c2->second;
}

void Face2D::Set1DRegion( IntField & ctrlpoints )
{
    this->ctrlpoints.push_back( ctrlpoints[ 0 ] );
    this->ctrlpoints.push_back( ctrlpoints[ 0 ] );
    this->ctrlpoints.push_back( ctrlpoints[ 1 ] );
    this->ctrlpoints.push_back( ctrlpoints[ 1 ] );
}

DomDataBasic::DomDataBasic()
{
    ;
}

DomDataBasic::~DomDataBasic()
{
}

DomData::DomData()
{
    ;
}

DomData::~DomData()
{
}

void DomData::ConstructCtrlPoint()
{
    std::map< int, IntSet >::iterator iter;
    for ( iter = this->pointToDomainMap.begin(); iter != this->pointToDomainMap.end(); ++ iter )
    {
        if ( iter->second.size() == 1 )
        {
            this->candidate_ctrlpoints.push_back( iter->first );
        }
    }

    int kkk = 1;
}

IntField & DomData::GetLinePoints( int line_id )
{
    int id = line_id - 1;
    return blkFaceSolver.lineList[ id ];
}

void DomData::ConstructBcPoint()
{
    IntField lineList;
    std::map< int, IntSet >::iterator iter;
    for ( iter = this->lineToDomainMap.begin(); iter != this->lineToDomainMap.end(); ++ iter )
    {
        if ( iter->second.size() == 1 )
        {
            lineList.push_back( iter->first );
        }
    }

    IntSet plist;

    for ( int iLine = 0; iLine < lineList.size(); ++ iLine )
    {
        int line_id = lineList[ iLine ];
        IntField & pointIdList = GetLinePoints( line_id );
        for ( int iPoint = 0; iPoint < pointIdList.size(); ++ iPoint )
        {
            int & pointId = pointIdList[ iPoint ];
            plist.insert( pointId );
        }
    }

    for ( IntSet::iterator iter = plist.begin(); iter != plist.end(); ++ iter )
    {
        this->candidate_bcpoints.push_back( *iter );
    }

    int kkk = 1;
}

void DomData::CalcDimBasic( int closedCurve )
{
    FindBcPointList2D( bcpointList );
    int nSeg = bcpointList.size();
    int nSegLoop = nSeg - 1;
    if ( closedCurve == 1 ) nSegLoop = nSeg;

    for ( int i = 0; i < nSegLoop; ++ i )
    {
        int i1 = i;
        int i2 = ( i + 1 ) % nSeg;
        int p1 = bcpointList[ i1 ];
        int p2 = bcpointList[ i2 ];
        IntField line( 2 );
        line[ 0 ] = p1;
        line[ 1 ] = p2;
        int line_id = blkFaceSolver.FindLineId( line );
        SegmentCtrl * segmentCtrl = line_Machine.GetSegmentCtrl( line_id + 1 );
        int dim = segmentCtrl->nPoint;
        this->bcdimList.push_back( dim );
    }

    int nCtrlPoint = this->ctrlpoints.size();

    IntField nij;
    int nCount = 0;
    int segCount = 0;
    int iCtrl = 0;
    for ( int i = 0; i < nSegLoop; ++ i )
    {
        int i1 = i;
        int i2 = ( i + 1 ) % nSeg;
        int p1 = bcpointList[ i1 ];
        int p2 = bcpointList[ i2 ];
        int dim = this->bcdimList[ i ];
        nCount += dim;
        segCount += 1;

        if ( IsCtrlPoint( p2 ) )
        {
            nCount -= ( segCount - 1 );
            nij.push_back( nCount );
            segCount = 0;
            nCount = 0;
            iCtrl ++;
        }
    }

    ni = nij[ 0 ];

    if ( nCtrlPoint > 2 )
    {
        nj = nij[ 1 ];
    }
    
}

void DomData::CalcDim2D()
{
    CalcDimBasic( 1 );
}

void DomData::CalcDim1D()
{
    CalcDimBasic( 0 );
}

void DomData::Normalize( int &d )
{
    if ( d != 0 )
    {
        d /= ABS( d );
    }
}

void DomData::FindBcPointList2D( IntField & bcpointList )
{
    int p1 = ctrlpoints[ 0 ];
    bcpointList.push_back( p1 );
    int prev = -1;
    int me = p1;
    int flag = -1;
    int next;

    while ( true )
    {
        FindNextPoint2D( bcpointList, prev, me, next, flag );
        if ( flag == 0 ) break;
        if ( next == p1 ) break;
        bcpointList.push_back( next );
        prev = me;
        me = next;
    };
    NormalBcPointList2D( bcpointList );
}

void DomData::NormalBcPointList2D( IntField & bcpointList )
{
    int p1 = ctrlpoints[ 0 ];
    int p2 = ctrlpoints[ 1 ];

    int nPoint = bcpointList.size();
    int dir = 1;
    for ( int i = 1; i < nPoint; ++ i )
    {
        int pt = bcpointList[ i ];
        bool flag = IsCtrlPoint( pt );
        if ( flag )
        {
            if ( p2 == pt )
            {
                dir = 1;
            }
            else
            {
                dir = -1;
            }
            break;
        }
    }

    if ( dir == -1 )
    {
        IntField tmpList = bcpointList;
        for ( int i = 1; i < nPoint; ++ i )
        {
            int j = nPoint - i;
            bcpointList[ i ] = tmpList[ j ];
        }
    }

}

bool DomData::IsCtrlPoint( int pt )
{
    int nPoint = this->ctrlpoints.size();
    for ( int i = 0; i < nPoint; ++ i )
    {
        int ip = this->ctrlpoints[ i ];
        if ( pt == ip )
        {
            return true;
        }
    }
    return false;
}

bool DomData::IsBcPoint( int pt )
{
    int nPoint = this->candidate_bcpoints.size();
    for ( int i = 0; i < nPoint; ++ i )
    {
        int ip = this->candidate_bcpoints[ i ];
        if ( pt == ip )
        {
            return true;
        }
    }
    return false;
}

void DomData::FindNextPoint2D( IntField & ptList, int prev, int me, int & next, int & flag )
{
    std::map< int, IntSet >::iterator iter;
    iter = this->pointToPointMap.find( me );
    IntSet & me_set = iter->second;
    flag = 0;
    for ( IntSet::iterator it = me_set.begin(); it != me_set.end(); ++ it )
    {
        next = * it;
        if ( IsBcPoint( next ) && ( next != prev ) && ( ! InArray( next, ptList ) ) )
        {
            flag = 1;
            break;
        }
    }
}

void DomData::CalcDomainCtrlPoints( IntField & blkControlpoints, IntField & localpt )
{
    for ( int i = 0; i < localpt.size(); ++ i )
    {
        int id = localpt[ i ] - 1;
        int pt = blkControlpoints[ id ];
        this->ctrlpoints.push_back( pt );
    }
}

void DomData::CalcDomainCtrlPoints( IntField & blk_ctrl_points )
{
    for ( int i = 0; i < blk_ctrl_points.size(); ++ i )
    {
        int pt = blk_ctrl_points[ i ];
        if ( InArray( pt, this->candidate_ctrlpoints ) )
        {
            this->ctrlpoints.push_back( pt );
        }
    }
}

bool DomData::IsBcLine( int line_id )
{
    std::map< int, IntSet >::iterator iter;
    iter = lineToDomainMap.find( line_id );
    return iter->second.size() == 1;
}

bool DomData::IsBcLine( IntSet &bclines, int line_id )
{
    IntSet::iterator iter;
    iter = bclines.find( line_id );
    return iter != bclines.end();
}

void DomData::RemoveBcLineId( IntSet &bclines, int line_id )
{
    bclines.erase( line_id );
}

void DomData::FindAllBoundaryLine( IntSet &bclines )
{
    std::map< int, IntSet >::iterator iter;
    for ( iter = lineToDomainMap.begin(); iter != lineToDomainMap.end(); ++ iter )
    {
        if ( iter->second.size() == 1 )
        {
            bclines.insert( iter->first );
        }
    }
}

bool DomData::FindNextBcPoint( int ps, int pt, int & pnext, IntSet &bclines )
{
    std::map< int, IntSet >::iterator iter;
    iter = this->pointToLineMap.find( pt );
    IntField lines;
    ONEFLOW::Set2Array( iter->second, lines );

    bool findflag = false;

    for ( int i = 0; i < lines.size(); ++ i )
    {
        int line_id = lines[ i ];
        if ( IsBcLine( bclines, line_id ) )
        {
            IntField pointIdList = GlobalGetLine( line_id );
            int p1 = pointIdList[ 0 ];
            int p2 = pointIdList[ 1 ];
            if ( p1 == pt )
            {
                pnext = p2;
            }
            else
            {
                pnext = p1;
            }
            RemoveBcLineId( bclines, line_id );
            findflag = true;
            break;
        }
    }
    return findflag && ( pnext != ps );
}

bool DomData::IsCornerPoints( int pt )
{
    for ( int i = 0; i < candidate_ctrlpoints.size(); ++ i )
    {
        int pp = candidate_ctrlpoints[ i ];
        if ( pp == pt )
        {
            return true;
        }
    }
    return false;
}

void DomData::CalcDomainCtrlPoints()
{
    //Start at any corner point
    //Using the point-to-edge data structure, find the edge.
    //If the edge is within the boundary (that is, the edge belongs to the boundary line), 
    //then move along the edge, find another point of the edge,
    //and remove the processed edge from the boundary. You can avoid going back.

    int ps = candidate_ctrlpoints[ 0 ];

    IntSet bclines;
    FindAllBoundaryLine( bclines );

    IntField bcpoints;

    int pt = ps;
    while ( true )
    {
        bcpoints.push_back( pt );
        int pnext = -1;
        bool flag = FindNextBcPoint( ps, pt, pnext, bclines );
        if ( ! flag ) break;
        pt = pnext;
    };

    for ( int i = 0; i < bcpoints.size(); ++ i )
    {
        int pct = bcpoints[ i ];
        if ( IsCornerPoints( pct ) )
        {
            this->ctrlpoints.push_back( pct );
        }
    }
}

void DomData::CalcBcCoor( CoorMap * coorMap, int closedCurve )
{
    int nSeg = bcpointList.size();
    int nSegLoop = nSeg - 1;
    if ( closedCurve == 1 ) nSegLoop = nSeg;

    int nCtrlPoint = this->ctrlpoints.size();
    int iCtrl = 0;
    for ( int i = 0; i < nSegLoop; ++ i )
    {
        int i1 = i;
        int i2 = ( i + 1 ) % nSeg;
        int p1 = bcpointList[ i1 ];
        int p2 = bcpointList[ i2 ];
        int dim = this->bcdimList[ i ];

        int cp1 = ctrlpoints[ iCtrl ];
        int cp2 = ctrlpoints[ ( iCtrl + 1 ) % nCtrlPoint ];

        CoorMap::iterator iter1 = coorMap->find( cp1 );
        CoorMap::iterator iter2 = coorMap->find( cp2 );
        int di = iter2->second.i - iter1->second.i;
        int dj = iter2->second.j - iter1->second.j;
        int dk = iter2->second.k - iter1->second.k;
        Normalize( di );
        Normalize( dj );
        Normalize( dk );

        int mi = iter1->second.i + di * ( dim - 1 );
        int mj = iter1->second.j + dj * ( dim - 1 );
        int mk = iter1->second.k + dk * ( dim - 1 );

        CalcCoor coor;
        coor.SetCoor( mi, mj, mk );
        coorMap->insert( std::pair<int, CalcCoor>( p2, coor ) );

        if ( IsCtrlPoint( p2 ) )
        {
            iCtrl ++;
        }
    }

    int kkk = 1;

}

void ConstructLineToDomainMap( int tid, IntField & idList, std::map< int, IntSet > & dataMap )
{
    for ( int i = 0; i < idList.size(); ++ i )
    {
        int sid = idList[ i ];
        ConstructInt2Map( sid, tid, dataMap );
    }
}

void ConstructIntList2Map( int tid, IntField & idList, std::map< int, IntSet > & dataMap )
{
    for ( int i = 0; i < idList.size(); ++ i )
    {
        int sid = idList[ i ];
        ConstructInt2Map( sid, tid, dataMap );
    }
}

void ConstructInt2Map( int sid, int tid, std::map< int, IntSet > & dataMap )
{
    std::map< int, IntSet >::iterator iter = dataMap.find( sid );

    if ( iter == dataMap.end() )
    {
        IntSet dataSet;
        dataSet.insert( tid );
        dataMap[ sid ] = dataSet;
    }
    else
    {
        iter->second.insert( tid );
    }
}

void ConstructPointToDomainMap( int tid, LinkField & pointIdLink, std::map< int, IntSet > & dataMap )
{
    int nLine = pointIdLink.size();

    for ( int iLine = 0; iLine < nLine; ++ iLine )
    {
        IntField & pointIdList = pointIdLink[ iLine ];
        ConstructIntList2Map( tid, pointIdList, dataMap );
    }
}

void ConstructPointToDomainMap( int tid, IntField & lineList, std::map< int, IntSet > & dataMap )
{
    for ( int iLine = 0; iLine < lineList.size(); ++ iLine )
    {
        int line_id = lineList[ iLine ] - 1;
        IntField & pointIdList = blkFaceSolver.lineList[ line_id ];

        ConstructIntList2Map( tid, pointIdList, dataMap );
    }
}

void ConstructPointToPointMap( LinkField & pointIdLink, std::map< int, IntSet > & dataMap )
{
    int nLine = pointIdLink.size();
    for ( int iLine = 0; iLine < nLine; ++ iLine )
    {
        IntField & pointIdList = pointIdLink[ iLine ];

        int & p1 = pointIdList[ 0 ];
        int & p2 = pointIdList[ 1 ];
        ConstructInt2Map( p1, p2, dataMap );
        ConstructInt2Map( p2, p1, dataMap );
    }
}

void ConstructPointToPointMap( IntField & lineList, std::map< int, IntSet > & dataMap )
{
    for ( int iLine = 0; iLine < lineList.size(); ++ iLine )
    {
        int line_id = lineList[ iLine ] - 1;
        IntField & pointIdList = blkFaceSolver.lineList[ line_id ];

        int & p1 = pointIdList[ 0 ];
        int & p2 = pointIdList[ 1 ];
        ConstructInt2Map( p1, p2, dataMap );
        ConstructInt2Map( p2, p1, dataMap );
    }
}

bool InArray( int ip, IntField & var_array )
{
    int nSize = var_array.size();
    for ( int i = 0; i < nSize; ++ i )
    {
        int pt = var_array[ i ];
        if ( pt == ip ) return true;
    }
    return false;
}

void GetPointIdLink( IntField & lineList, LinkField & pointIdLink )
{
    for ( int iLine = 0; iLine < lineList.size(); ++ iLine )
    {
        int line_id = lineList[ iLine ] - 1;
        IntField & pointIdList = blkFaceSolver.lineList[ line_id ];
        pointIdLink.push_back( pointIdList );
    }
}

void GetUnitInt( int &d )
{
    if ( d != 0 )
    {
        d /= ABS( d );
    }
}

void GetUnitDir( CalcCoor & c )
{
    GetUnitInt( c.i );
    GetUnitInt( c.j );
    GetUnitInt( c.k );
}

EndNameSpace

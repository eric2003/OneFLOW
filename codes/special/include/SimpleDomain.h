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


#pragma once
#include "HXDefine.h"
#include "Mid.h"
#include "CalcCoor.h"
#include <set>
#include <map>
#include <fstream>
using namespace std;
BeginNameSpace( ONEFLOW )

class Block3D;
class SDomain;
class CalcCoor;
class Face2D;

class BlkF2C
{
public:
    BlkF2C();
    ~BlkF2C();
public:
    int id, type, bctype;
    IntField cellList;
    IntField posList;
};

class Face2D
{
public:
    Face2D();
    ~Face2D();
public:
    int face_id;
    int bcType;
    IntField ctrlpoints;
    CalcCoor st, ed;
    CalcCoor p1, p2;
    Face2D * t;
public:
    void CalcRegion();
    void CalcStEd( CoorMap * coorMap );
    void Set1DRegion( IntField & ctrlpoints );
};

class DomDataBasic
{
public:
    DomDataBasic();
    ~DomDataBasic();
public:
    IntField candidate_ctrlpoints;
    IntField candidate_bcpoints;
    IntField ctrlpoints;
    IntField original_ctrlpoints;
    IntField bcpointList;
    IntField bcdimList;

    int ni, nj;
public:
    map< int, IntSet > pointToDomainMap;
    map< int, IntSet > pointToPointMap;
    map< int, IntSet > lineToDomainMap;
};

class DomData : public DomDataBasic
{
public:
    DomData();
    ~DomData();
public:
    IntField & GetLinePoints( int line_id );
    void ConstructCtrlPoint();
    void ConstructBcPoint();
    void CalcDimBasic( int closedCurve );
    void CalcDim2D();
    void CalcDim1D();
    void CalcBcCoor( CoorMap * coorMap, int iloop );
    void Normalize( int &d );
    void FindBcPointList2D( IntField & bcpointList );
    void NormalBcPointList2D( IntField & bcpointList );
    void FindNextPoint2D( IntField & ptList, int prev, int me, int & next, int & flag );
    bool IsBcPoint( int pt );
    bool IsCtrlPoint( int pt );
    void CalcDomainCtrlPoints( IntField & blkControlpoints, IntField & localpt );
    void CalcDomainCtrlPoints( IntField & blk_ctrl_points );
};


void ConstructInt2Map( int sid, int tid, map< int, IntSet > & dataMap );
void ConstructIntList2Map( int tid, IntField & idList, map< int, IntSet > & dataMap );
void ConstructLineToDomainMap( int tid, IntField & idList, map< int, IntSet > & dataMap );
void ConstructPointToDomainMap( int tid, IntField & lineList, map< int, IntSet > & dataMap );
void ConstructPointToPointMap( IntField & lineList, map< int, IntSet > & dataMap );

void ConstructPointToDomainMap( int tid, LinkField & pointIdLink, map< int, IntSet > & dataMap );
void ConstructPointToPointMap( LinkField & pointIdLink, map< int, IntSet > & dataMap );
bool InArray( int ip, IntField & var_array );

void GetPointIdLink( IntField & lineList, LinkField & pointIdLink );

void GetUnitInt( int &d );
void GetUnitDir( CalcCoor & c );

EndNameSpace
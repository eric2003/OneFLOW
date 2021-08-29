/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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
#include "PointFactory.h"
#include <map>
using namespace std;

BeginNameSpace( ONEFLOW )

class PointAction
{
public:
    PointAction();
    ~PointAction();
public:
    typedef Point< Real > PointType;
    typedef map< PointType , int, PointCompare< Real > > PointMap;
public:
    PointMap pointMap;
    HXVector< PointType > pointList;
public:
    UInt GetNPoint() { return pointList.size(); }
    int AddPoint( Real xm, Real ym, Real zm );
    int DeletePoint( Real xm, Real ym, Real zm );
    int DeletePoint( PointAction::PointType & point );
    int FindPoint( Real xm, Real ym, Real zm );
    void ModifyPointIndexAfterDelete( int pid );
protected:
    int AddPoint( PointAction::PointType & point );
    int FindPointId( PointAction::PointType & point );
    bool FindPoint( PointAction::PointType & point, PointAction::PointMap::iterator & iter );
};

class Ids
{
public:
    Ids();
    ~Ids();
public:
    vector< int > ids;
    vector< int > sorted_ids;
    int type;
};

class CompareIds
{
public:
    bool operator()( const Ids & lhs, const Ids & rhs ) const;
};

class IdTool
{
public:
    IdTool() ;
    ~IdTool();
public:
   typedef map< Ids, int, CompareIds > IDSMap;
   IDSMap ids_map;
   vector< Ids > ids_list;
   Ids vint;
public:
    bool NotFind( IDSMap::iterator & iter );
    IDSMap::iterator FindIds( const vector< int > & ids, int type );
    int AddIds( vector< int > & ids, int type );
    int AddData();
    void ModifyDataIndex( const Ids &var, int new_id );
};

class UnitElement;

class TopoSort
{
public:
    TopoSort();
    ~TopoSort();
public:
    IdTool faceIdTool;
    vector< int > real_face;
public:
    vector< int > lc, rc;
    vector< int > lc_pos, rc_pos;
    vector< int > fTypes;
    vector< int > fBcTypes;
    vector< int > bcTypes;
    int nCells;
public:
    void GetElementFace( UnitElement * unitElement, vector< int > & element, int facePos, vector< int > & face, int & faceType );
    void AddSingleFace( UnitElement * unitElement, vector< int > & element, int facePos, int  iCell );
    void ModifyFace( int face_id, int iCell, int face_pos );
    void AddNewFace( int iCell, int face_pos, int faceType );
    void AddElementFaces( vector< int > & element, int eType, int iCell );
public:
    void ReorderFaces();
    void CalcOrderMap( vector<int > & orderMap );
    void ReOrder( vector< int > & varList, vector< int > & orderMap );
    void ReOrderMapdata( IdTool & faceIdTool, vector< int > & orderMap );
    void TopoPostprocess();
public:
    void ScanBcFace();
    void SetBcGhostCell();
};

class TopoAction
{
public:
    TopoAction();
    ~TopoAction();
public:
    void AddElement( int p1, int p2, int eType );
public:
    void CalcTopology();
public:
    typedef map< Ids, int, CompareIds > PointMap;
public:
    IdTool elementIdTool;
public:
    TopoSort * topo_sort;
};

class MyBcRegion
{
public:
    MyBcRegion();
    ~MyBcRegion();
public:
    string name;

};


class BcAction
{
public:
    BcAction();
    ~BcAction();
public:
    void CreatBCRegion( const string "LeftOutFlow", ONEFLOW::BCOutflow );
};


class SmartGrid
{
public:
    SmartGrid();
    ~SmartGrid();
public:
    PointAction * point_action;
    TopoAction * topo_action;
    BcAction * topo_action;
public:
    int AddPoint( Real x, Real y, Real z );
    void TestAddDeletePoints();
    void Run();
public:
    void AddElement( int p1, int p2, int eType );
    void GenerateGrid( int ni, Real xmin, Real xmax );
    void CalcTopology();
    void TopoPostprocess();
    void ReorderFaces();
    void CalcOrderMap( vector< int > & orderMap );
};


EndNameSpace
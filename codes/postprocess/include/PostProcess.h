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
#include "Configure.h"
#include "HXType.h"
#include "HXDefine.h"
#include "HXSort.h"
#include <vector>
#include <set>
using namespace std;

BeginNameSpace( ONEFLOW )

class MakeCurveClass;

class Post
{
public:
    Post();
    ~Post();
public:
    void Run();
};

class VectDir
{
public:
    VectDir();
    ~VectDir();
public:
    Real dx, dy;
public:
    void Normalize();
    Real Angle( VectDir * vd );
};

void GetCurveTitle( StringField & title );
void GetCurveValueTitle( StringField & title );

class PointEdgeClass
{
public:
    PointEdgeClass();
    ~PointEdgeClass();
public:
    set<int> idset;
    vector<int> idlist;
public:
    void EdgeToPoint( vector< IntField > & edgeList );
};

class CurveData
{
public:
    CurveData();
    ~CurveData();
public:
    vector< Real > xList, yList, zList;
    vector< Real > cpList, cfList;
    vector< int > extremeList;
    vector< int > idList;
    vector< HXSort< Real > > xsList;
    vector< Real > dsList;
    vector< VectDir > vectDirList;
    set< int > idset;
    vector< int > newidlist;
    vector< set< int > > p2p;
    vector< IntField > extremeEdgeList;
    set < HXSort< IntField > > searchEdgeList;
    PointEdgeClass pec;
    string file_prestr;
public:
    void SetFilePreStr( const string & file_prestr );
    void AddPoint( Real x, Real y, Real z );
    void AddValue( Real cp, Real cf );
    void AddPointValue( Real x, Real y, Real z, Real cp, Real cf );
    void FindExtremePoint();
    bool InTriangle( int p, int q, int r, int s );
    bool ToLeft( int p, int q, int s );
    void FindExtremeEdge();
    bool CheckEdge( int p, int q );
    void AddExtremeEdge( int p, int q );
    bool FindEdge( int p, int q );
    void InitP2p();
public:
    void AddVector( int p );
    void RemovePoint( int p );
    int FindNextExtremeEdgePoint( int p1, int p2 );
    bool FindNearestPoint( int p0, int &pNearest );
    int FindXminIndex();
    void CreateXSortList();
    void SetExtremePoint();
    void ConstructFinalList();
    IntField GetFistEdge();
    bool FindNextPoint( int p1, int & p2, VectDir * vd );
public:
    void VisualCreateList();
    void VisualCreateList( PointEdgeClass * pec );
    void VisualCreateList( vector<int> & newidlist, const string& fileName );
    void VisualCurveValue();
    void VisualCurveValue( vector<int> & newidlist, const string& fileName );
    void Visual( const string & fileName, vector< Real > & x, vector< Real > & y, vector< Real > & z );
    void VisualCurveValue( const string & fileName, vector< Real > & x, vector< Real > & y, vector< Real > & z, vector< Real > & cp, vector< Real > & cf );
    void Visual( const string & fileName );

};

class MakeCurveClass
{
public:
    MakeCurveClass();
    ~MakeCurveClass();
public:
    CurveData totalPart, mainPart, slapPart;
    CurveData upperMain, lowerMain;
public:
    void Run();
public:
    void SplitCurve();
    void NLR7301MainSplit();
    void AddPoint( Real x, Real y, Real z );
    void AddPointValue( Real x, Real y, Real z, Real cp, Real cf );
    void Visual();
};


void PostSimu();

EndNameSpace

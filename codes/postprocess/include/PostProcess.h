/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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
    std::set<int> idset;
    std::vector<int> idlist;
public:
    void EdgeToPoint( std::vector< IntField > & edgeList );
};

class CurveData
{
public:
    CurveData();
    ~CurveData();
public:
    std::vector< Real > xList, yList, zList;
    std::vector< Real > cpList, cfList;
    std::vector< int > extremeList;
    std::vector< int > idList;
    std::vector< HXSort< Real > > xsList;
    std::vector< Real > dsList;
    std::vector< VectDir > vectDirList;
    std::set< int > idset;
    std::vector< int > newidlist;
    std::vector< std::set< int > > p2p;
    std::vector< IntField > extremeEdgeList;
    std::set< HXSort< IntField > > searchEdgeList;
    PointEdgeClass pec;
    std::string file_prestr;
public:
    void SetFilePreStr( const std::string & file_prestr );
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
    void VisualCreateList( std::vector<int> & newidlist, const std::string& fileName );
    void VisualCurveValue();
    void VisualCurveValue( std::vector<int> & newidlist, const std::string& fileName );
    void Visual( const std::string & fileName, std::vector< Real > & x, std::vector< Real > & y, std::vector< Real > & z );
    void VisualCurveValue( const std::string & fileName, std::vector< Real > & x, std::vector< Real > & y, std::vector< Real > & z, std::vector< Real > & cp, std::vector< Real > & cf );
    void Visual( const std::string & fileName );

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

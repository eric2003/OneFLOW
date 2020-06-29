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
#include "CalcCoor.h"
#include "SimpleDomain.h"
#include <map>
using namespace std;

BeginNameSpace( ONEFLOW )

class SDomain;

class SLine
{
public:
    SLine();
    ~SLine();
public:
    int line_id;
    int ni;
    RealField x1d, y1d, z1d;
    IntField ctrlpoints;
public:
    void SetDomainBcMesh( SDomain * sDomain );
    void ConstructCtrlPoints();
    void Alloc();
    void CopyMesh();

};

class MLine : public DomData
{
public:
    MLine();
    ~MLine();
public:
    int pos;
    IntField lineList;
    HXVector< SLine * > slineList;
    CoorMap * coorMap;
public:
    map< int, IntSet > pointToLine;
public:
    void ConstructLineToDomainMap();
    void ConstructLineToDomainMap( int domain_id, map< int, IntSet > & lineToDomainMap );
    void ConstructPointToDomainMap();
    void ConstructPointToDomainMap( int domain_id, map< int, IntSet > & pointToDomainMap );
    void ConstructPointToPointMap();
    void ConstructPointToPointMap( map< int, IntSet > & pointToPointMap );
public:
    void AddSubLine( int line_id );
    void ConstructDomainTopo();
    void ConstructCtrlPoint();
    void ConstructSLineCtrlPoint();
    void CalcCoor( CoorMap * localCoorMap );
    void SetDomainBcMesh( SDomain * sDomain );
    void CreateInpFaceList( HXVector< Face2D * > &facelist );
};


EndNameSpace
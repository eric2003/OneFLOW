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
#include "HXCgns.h"

BeginNameSpace( ONEFLOW )

#ifdef ENABLE_CGNS

class CgnsZone;
class NodeMesh;

class CgnsCoor
{
public:
    CgnsCoor( CgnsZone * cgnsZone );
    ~CgnsCoor();
public:
    int ndim;
    int nCoor;
    IntField nNodeList;
    HXVector< DataType_t > typeList;
    HXVector< void * > coor;
    CgnsZone * cgnsZone;
    NodeMesh * nodeMesh;
public:
    CgInt irmin[ 3 ], irmax[ 3 ], cellSize[ 3 ];
protected:
    CgInt nNode, nCell;
public:
    CgInt GetNNode();
    CgInt GetNCell();

    void SetNNode( CgInt nNode );
    void SetNCell( CgInt nCell );
public:
    void * GetCoor( int iCoor ) { return coor[ iCoor ]; };
    void SetAllData( RealField & x, RealField & y, RealField & z );
    void Alloc( int iCoor, int nNode, DataType_t data_type );
public:
    void SetData( int iCoor, DataType_t data_type, Real * var );
    void DeAlloc();
public:
    void ReadCgnsGridCoordinates();
    void FreeMesh();
public:
    NodeMesh * GetNodeMesh();
    void SetDimension();
    void SetDimension( CgnsCoor * cgnsCoorIn );
    void SetDimensionStr();
    void SetDimensionUns();
};


#endif

EndNameSpace
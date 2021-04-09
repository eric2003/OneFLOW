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
#include "ScalarGrid.h"
#include "metis.h"
#include "ScalarIFace.h"
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

typedef vector<idx_t> MetisIntList;
class ScalarGrid;

class MetisPart
{
public:
    MetisPart( ScalarGrid * ggrid );
    ~MetisPart();
public:
    MetisIntList xadj;
    MetisIntList adjncy;
    ScalarGrid * ggrid;
public:
    void MetisPartition( int nPart, MetisIntList & cellzone );
private:
    void ScalarGetXadjAdjncy( ScalarGrid * ggrid, MetisIntList & xadj, MetisIntList & adjncy );
    void ScalarPartitionByMetis( idx_t nCells, MetisIntList & xadj, MetisIntList & adjncy, int nPart, MetisIntList & cellzone );

};

class NetGrid
{
public:
    NetGrid();
    ~NetGrid();
public:
    vector< ScalarGrid * > grids;
public:
    void AllocateGrid( int nZones );
    void DeAllocateGrid();
};

class GridTopo
{
public:
    GridTopo();
    ~GridTopo();
public:
    vector<int> faceid;
    vector<int> nodeid;
    vector<int> facetype;
    set<int> nodeset;
    IntList lc;
    IntList rc;
    int zoneid;
public:
    void AddFaceId( int iFace );
    void AddFaceType( int faceType );
    void ReconstructNode( EList & faces );
};

class GridTopos
{
public:
    GridTopos();
    ~GridTopos();
public:
    vector< GridTopo > data;
    vector<int> gLCells;
    ScalarIFaces scalarIFaces;
public:
    void Allocate( int nZones );
    GridTopo & operator [] ( int i )
    {
        return data[ i ];
    }
public:
    void CalcGlobal2LocalCells( MetisIntList & cellzone );
    void CalcInterface();
};

class Part
{
public:
    Part();
    ~Part();
public:
    ScalarGrid * ggrid;
    NetGrid * netGrid;
    MetisIntList cellzone;
    int nPart;
    GridTopos gtopos;
public:
    void PartitionGrid( ScalarGrid * ggrid, int nPart, NetGrid * netGrid );
    void CalcCellZone();
    void ReconstructAllZones();
    void CalcGlobalInterface();
    void ReconstructGrid();
    void ReconstructNode();
    
};


EndNameSpace
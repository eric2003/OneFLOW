/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
	Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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
#include "GridDef.h"
#include "HXCgns.h"
#include <vector>
#include <string>
#include <fstream>
using namespace std;

BeginNameSpace( ONEFLOW )

class StrGrid;
class NodeMesh;
class Grid;
class CgnsMultiBase;
class CgnsZone;
class GridElem;
class Su2Grid;

#ifdef ENABLE_CGNS

class CgnsFactory
{
public:
	CgnsFactory();
	~CgnsFactory();
public:
    CgnsMultiBase * cgnsMultiBase;

    Grids cmpGrids;
    HXVector< GridElem * > gridElems;

    int nZone;
    int nOriZone;
public:
    void GenerateGrid();
    void ReadCgnsGrid();
    void ProcessGrid();
public:
    void CommonToOneFlowGrid();
    void CommonToUnsGrid();
    void CommonToStrGrid();
public:
    CgnsZone * GetCreateZone( int cgnsZoneId = 0 );
    void MergeToSingleZone( Grids & grids, HXVector< Int3D * > & unsIdList, NodeMesh * nodeMesh, int & nNode, int & nCell );
    void PrepareSection( Grids & grids, int cgnsZoneId );
    void PrepareSectionBasic( Grids & grids, CgnsZone * cgnsZone );
    void FillSection( Grids & grids, HXVector< Int3D * > & unsIdList );
public:
    void PrepareSU2Section( Su2Grid* su2Grid, int cgnsZoneId );
    void PrepareSU2SectionBasic( Su2Grid* su2Grid, CgnsZone * cgnsZone );
    void MergeSU2ToSingleZone( Su2Grid* su2Grid, HXVector< Int3D * > & unsIdList, NodeMesh * nodeMesh, int & nNode, int & nCell );
    void FillSU2Section( Su2Grid* su2Grid );
public:
    void CgnsStr2Uns( Grid *& grid, int zId );
    void ConvertStrCgns2UnsCgnsGrid();
    void AllocateGridElem();
    void DeAllocateGridElem();
    void PrepareUnsCompGrid();
    void AllocateCmpGrid();

	//转换为oneflow计算所用的网格
    void GenerateCmpGrid();

protected:
    void GenerateStrCmpGrid();
    void GenerateUnsCmpGrid();
};

class Grid;
class PointSearch;
class NodeMesh;
void ComputeUnsId( StrGrid * grid, PointSearch * pointSearch, Int3D * unsId );

int OneFlow2CgnsZoneType( int zoneType );
int Cgns2OneFlowZoneType( int zoneType );

class BcRegion;
void SetUnsBcConn( BcRegion * bcRegion, CgIntField& conn, int & pos, Int3D & unsId );

#endif

EndNameSpace
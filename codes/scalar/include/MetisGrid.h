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
#include "ScalarGrid.h"
#include "metis.h"
#include <vector>
#include <set>
#include <map>


BeginNameSpace( ONEFLOW )

typedef std::vector<idx_t> MetisIntList;
class ScalarGrid;

class MetisSplit
{
public:
    MetisSplit();
    ~MetisSplit();
public:
    MetisIntList xadj;
    MetisIntList adjncy;
public:
    void MetisPartition( ScalarGrid * ggrid, int nPart, MetisIntList & cellzone );
    void ManualPartition( ScalarGrid * ggrid, int nPart, MetisIntList & cellzone );
private:
    void ScalarGetXadjAdjncy( ScalarGrid * ggrid, MetisIntList & xadj, MetisIntList & adjncy );
    void ScalarPartitionByMetis( idx_t nCells, MetisIntList & xadj, MetisIntList & adjncy, int nPart, MetisIntList & cellzone );

};

class ScalarIFace;

class GridPartition
{
public:
    GridPartition();
    ~GridPartition();
public:
    ScalarGrid * ggrid;
    int nPart;
    std::vector< ScalarGrid * > * grids;
public:
    int GetNZones();
    void AllocateGrid( int nZones );
    void PartitionGrid( ScalarGrid * ggrid, int nPart, std::vector< ScalarGrid * > * grids );
    void ReconstructGridFaceTopo();
    void ReconstructInterfaceTopo();
    void ReconstructNode();
    void ReconstructNeighbor();
    void CalcInterfaceToBcFace();
};


EndNameSpace

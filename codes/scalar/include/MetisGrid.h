/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
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
#include "NamespaceMacros.h"
#include "HXType.h"
#include "HXDefine.h"
#include "ScalarGrid.h"
#include "metis.h"
#include <memory>
#include <vector>
#include <set>
#include <map>


BeginNameSpace( ONEFLOW )

using MetisIntList = std::vector<idx_t>;
class ScalarGrid;

class MetisSplit
{
public:
    MetisIntList MetisPartition( const ScalarGrid & ggrid, int nPart );
    MetisIntList ManualPartition( const ScalarGrid & ggrid );
private:
    std::pair< MetisIntList, MetisIntList > ScalarGetXadjAdjncy( const ScalarGrid & ggrid );
    MetisIntList ScalarPartitionByMetis( idx_t nCells, const MetisIntList & xadj, const MetisIntList & adjncy, int nPart );

};

class ScalarIFace;

class GridPartition
{
public:
    std::vector< std::unique_ptr< ScalarGrid > > PartitionGrid( const ScalarGrid & ggrid, int nPart );
private:
    std::vector< std::unique_ptr< ScalarGrid > > AllocateGrid( int nZones );
    std::vector< std::unique_ptr< ScalarGrid > > ReconstructGridFaceTopo( const ScalarGrid & ggrid, int nPart );
    void ReconstructInterfaceTopo( std::vector< std::unique_ptr< ScalarGrid > > & grids );
    void ReconstructNode( const ScalarGrid & ggrid, std::vector< std::unique_ptr< ScalarGrid > > & grids );
    void ReconstructNeighbor( std::vector< std::unique_ptr< ScalarGrid > > & grids );
    void CalcInterfaceToBcFace( std::vector< std::unique_ptr< ScalarGrid > > & grids );
};


EndNameSpace

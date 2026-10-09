#include <memory>
#include <utility>
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

#include "MetisGrid.h"
#include "ScalarGrid.h"
#include "ScalarIFace.h"
#include "Constant.h"
#include "HXCgns.h"
#include "ElementHome.h"
#include "HXMath.h"
#include "Boundary.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <limits>


BeginNameSpace( ONEFLOW )

MetisIntList MetisSplit::ManualPartition( const ScalarGrid & ggrid )
{
	int nCells = ggrid.GetNCells();

	MetisIntList cellzone( nCells );
	std::vector< int > tmp;
	for ( int iCell = 0; iCell < nCells; iCell += 2 )
	{
		tmp.push_back( iCell );
	}

	for ( int iCell = 1; iCell < nCells; iCell += 2 )
	{
		tmp.push_back( iCell );
	}

	for ( int iCell = 0; iCell < nCells; ++ iCell )
	{
		cellzone[ iCell ] = tmp[ iCell ];
	}

	return cellzone;
}

MetisIntList MetisSplit::MetisPartition( const ScalarGrid & ggrid, int nPart )
{
	const int nCells = ggrid.GetNCells();
	if ( nCells <= 0 )
	{
		throw std::invalid_argument( "MetisSplit::MetisPartition: input grid must contain at least one cell" );
	}
	if ( nPart <= 0 || nPart > nCells )
	{
		throw std::invalid_argument( "MetisSplit::MetisPartition: nPart must be between 1 and the number of cells" );
	}

	if ( nPart == nCells )
	{
		return ManualPartition( ggrid );
	}

	auto graph = ScalarGetXadjAdjncy( ggrid );
	return ScalarPartitionByMetis( nCells, graph.first, graph.second, nPart );
}

std::pair< MetisIntList, MetisIntList > MetisSplit::ScalarGetXadjAdjncy( const ScalarGrid & ggrid )
{
	const int nCells = ggrid.GetNCells();
	const int nFaces = ggrid.GetNFaces();
	const int nBFaces = ggrid.GetNBFaces();
	if ( nCells <= 0 || nBFaces < 0 || nFaces < nBFaces ||
		 nBFaces > std::numeric_limits< int >::max() - nCells )
	{
		throw std::invalid_argument( "MetisSplit::ScalarGetXadjAdjncy: invalid grid topology counts" );
	}
	if ( ggrid.lc.GetNElements() != static_cast< size_t >( nFaces ) ||
		 ggrid.rc.GetNElements() != static_cast< size_t >( nFaces ) ||
		 ggrid.bcTypes.GetNElements() != static_cast< size_t >( nBFaces ) )
	{
		throw std::runtime_error( "MetisSplit::ScalarGetXadjAdjncy: face topology arrays have inconsistent sizes" );
	}

	// Validate physical cell references before CalcC2C indexes the adjacency rows.
	for ( int iFace = 0; iFace < nBFaces; ++ iFace )
	{
		const int leftCell = ggrid.lc[ iFace ];
		if ( leftCell < 0 || leftCell >= nCells )
		{
			throw std::runtime_error( "MetisSplit::ScalarGetXadjAdjncy: boundary face references an invalid physical cell" );
		}
	}
	for ( int iFace = nBFaces; iFace < nFaces; ++ iFace )
	{
		const int leftCell = ggrid.lc[ iFace ];
		const int rightCell = ggrid.rc[ iFace ];
		if ( leftCell < 0 || leftCell >= nCells || rightCell < 0 || rightCell >= nCells )
		{
			throw std::runtime_error( "MetisSplit::ScalarGetXadjAdjncy: internal face references an invalid physical cell" );
		}
	}

	EList c2c;
	ggrid.CalcC2C( c2c );
	if ( c2c.GetNElements() != static_cast< size_t >( nCells ) )
	{
		throw std::runtime_error( "MetisSplit::ScalarGetXadjAdjncy: cell adjacency row count does not match the cell count" );
	}

	MetisIntList xadj( static_cast< size_t >( nCells ) + 1 );
	MetisIntList adjncy;
	const size_t nInternalFaces = static_cast< size_t >( nFaces - nBFaces );
	if ( nInternalFaces <= adjncy.max_size() / 2 )
	{
		adjncy.reserve( 2 * nInternalFaces );
	}

	xadj[ 0 ] = 0;
	for ( int iCell = 0; iCell < nCells; ++ iCell )
	{
		for ( const int neighbor : c2c[ iCell ] )
		{
			if ( neighbor < 0 || neighbor >= nCells + nBFaces )
			{
				throw std::runtime_error( "MetisSplit::ScalarGetXadjAdjncy: cell adjacency contains an invalid cell index" );
			}
			// METIS partitions physical cells only; interface ghost cells are not graph vertices.
			if ( neighbor >= nCells )
			{
				continue;
			}
			if ( adjncy.size() >= static_cast< size_t >( std::numeric_limits< idx_t >::max() ) )
			{
				throw std::overflow_error( "MetisSplit::ScalarGetXadjAdjncy: adjacency exceeds METIS index range" );
			}
			adjncy.push_back( static_cast< idx_t >( neighbor ) );
		}
		xadj[ iCell + 1 ] = static_cast< idx_t >( adjncy.size() );
	}

	return { std::move( xadj ), std::move( adjncy ) };
}

MetisIntList MetisSplit::ScalarPartitionByMetis( idx_t nCells, const MetisIntList & xadj, const MetisIntList & adjncy, int nPart )
{
	if ( nCells <= 0 )
	{
		throw std::invalid_argument( "MetisSplit::ScalarPartitionByMetis: number of cells must be positive" );
	}
	if ( nPart <= 0 || nPart > nCells )
	{
		throw std::invalid_argument( "MetisSplit::ScalarPartitionByMetis: nPart must be between 1 and the number of cells" );
	}
	if ( xadj.size() != static_cast< size_t >( nCells ) + 1 || xadj.empty() || xadj[ 0 ] != 0 )
	{
		throw std::invalid_argument( "MetisSplit::ScalarPartitionByMetis: invalid CSR row offsets" );
	}
	if ( adjncy.size() > static_cast< size_t >( std::numeric_limits< idx_t >::max() ) )
	{
		throw std::invalid_argument( "MetisSplit::ScalarPartitionByMetis: adjacency array exceeds METIS index range" );
	}
	for ( idx_t iCell = 0; iCell < nCells; ++ iCell )
	{
		if ( xadj[ iCell ] < 0 || xadj[ iCell + 1 ] < xadj[ iCell ] ||
			 static_cast< size_t >( xadj[ iCell + 1 ] ) > adjncy.size() )
		{
			throw std::invalid_argument( "MetisSplit::ScalarPartitionByMetis: CSR row offsets are not monotonic or exceed adjacency data" );
		}
	}
	if ( static_cast< size_t >( xadj[ nCells ] ) != adjncy.size() )
	{
		throw std::invalid_argument( "MetisSplit::ScalarPartitionByMetis: final CSR offset does not match adjacency size" );
	}
	for ( const idx_t neighbor : adjncy )
	{
		if ( neighbor < 0 || neighbor >= nCells )
		{
			throw std::invalid_argument( "MetisSplit::ScalarPartitionByMetis: adjacency contains an invalid cell index" );
		}
	}

	MetisIntList cellzone( static_cast< size_t >( nCells ) );
	MetisIntList metisXadj = xadj;
	MetisIntList metisAdjncy = adjncy;
	idx_t ncon = 1;
	idx_t * vwgt = 0;
	idx_t * vsize = 0;
	idx_t * adjwgt = 0;
	float * tpwgts = 0;
	float * ubvec = 0;
	idx_t options[ METIS_NOPTIONS ];
	idx_t wgtflag = 0;
	idx_t numflag = 0;
	idx_t objval;
	idx_t nZone = nPart;
	idx_t emptyAdjacency = 0;
	idx_t * xadjData = metisXadj.data();
	idx_t * adjncyData = metisAdjncy.empty() ? & emptyAdjacency : metisAdjncy.data();

	const int optionsStatus = METIS_SetDefaultOptions( options );
	if ( optionsStatus != METIS_OK )
	{
		throw std::runtime_error( "MetisSplit::ScalarPartitionByMetis: failed to initialize METIS options" );
	}

	std::cout << "Now begining partition graph!\n";
	int partitionStatus = METIS_OK;
	if ( nZone > 8 )
	{
		std::cout << "Using K-way Partitioning!\n";
		partitionStatus = METIS_PartGraphKway( & nCells, & ncon, xadjData, adjncyData, vwgt, vsize, adjwgt,
			& nZone, tpwgts, ubvec, options, & objval, cellzone.data() );
	}
	else
	{
		std::cout << "Using Recursive Partitioning!\n";
		partitionStatus = METIS_PartGraphRecursive( & nCells, & ncon, xadjData, adjncyData, vwgt, vsize, adjwgt,
			& nZone, tpwgts, ubvec, options, & objval, cellzone.data() );
	}

	if ( partitionStatus != METIS_OK )
	{
		throw std::runtime_error( "MetisSplit::ScalarPartitionByMetis: METIS graph partitioning failed" );
	}

	std::cout << "The interface number: " << objval << std::endl;
	std::cout << "Partition is finished!\n";
	return cellzone;
}

std::vector< std::unique_ptr< ScalarGrid > > GridPartition::PartitionGrid( const ScalarGrid & ggrid, int nPart )
{
    std::vector< std::unique_ptr< ScalarGrid > > grids;

    grids = ReconstructGridFaceTopo( ggrid, nPart );
    ReconstructNeighbor( grids );
    ReconstructInterfaceTopo( grids );
    CalcInterfaceToBcFace( grids );
    ReconstructNode( ggrid, grids );

    return grids;
}

std::vector< std::unique_ptr< ScalarGrid > > GridPartition::AllocateGrid( int nZones )
{
    std::vector< std::unique_ptr< ScalarGrid > > grids;
    for ( int iZone = 0; iZone < nZones; ++ iZone )
    {
        auto grid = std::make_unique< ScalarGrid >();
        grid->id = iZone;
        grids.push_back( std::move( grid ) );
    }
    return grids;
}

std::vector< std::unique_ptr< ScalarGrid > > GridPartition::ReconstructGridFaceTopo( const ScalarGrid & ggrid, int nPart )
{
	// Calculate the cell-to-zone mapping before allocating zone-local topology.
	MetisIntList cellzone = MetisSplit::MetisPartition( ggrid, nPart );
	const int nCells = ggrid.GetNCells();
	if ( cellzone.size() != static_cast< size_t >( nCells ) )
	{
		throw std::runtime_error( "GridPartition::ReconstructGridFaceTopo: partition result size does not match the cell count" );
	}
	for ( int iCell = 0; iCell < nCells; ++ iCell )
	{
		if ( cellzone[ iCell ] < 0 || cellzone[ iCell ] >= nPart )
		{
			throw std::runtime_error( "GridPartition::ReconstructGridFaceTopo: partition result contains an invalid zone id" );
		}
	}

	std::vector< std::unique_ptr< ScalarGrid > > grids = AllocateGrid( nPart );

	int nZones = static_cast< int >( grids.size() );
	int nFaces = ggrid.GetNFaces();
	int nBFaces = ggrid.GetNBFaces();

	std::vector<int> zoneCount( nZones, 0 );
	std::vector<int> localCells; //global cell id -> local cell id
	localCells.resize( nCells );

	for ( int iCell = 0; iCell < nCells; ++ iCell )
	{
		int iZone = cellzone[ iCell ];
		localCells[ iCell ] = zoneCount[ iZone ] ++;
		int eType = ggrid.eTypes[ iCell ];
		ScalarGrid & grid = *grids[ iZone ];
		grid.eTypes.AddData( eType );
	}

	//First scan the global physical boundary
	for ( int iFace = 0; iFace < nBFaces; ++ iFace )
	{
		int lc = ggrid.lc[ iFace ];
		int bctype = ggrid.bcTypes[ iFace ];
		int lZone = cellzone[ lc ];
		//global face id = iFace, local face id = faceid.size();
		//global face node: 20,10,30,40, local face node 1 2 4 3 for example
		//global coor x[20],y[20],z[20],x[10],y[10],z[10]
		//local coor x[1],y[1],z[1],x[2],y[2],z[2]
		int localCell = localCells[ lc ];
		int ftype = ggrid.fTypes[ iFace ];
		ScalarGrid & gridL = *grids[ lZone ];
		gridL.AddFaceType( ftype );
		gridL.AddPhysicalBcFace( iFace, bctype, localCell, ONEFLOW::INVALID_INDEX );
	}

	//Then scan the internal block interface
	for ( int iFace = nBFaces; iFace < nFaces; ++ iFace )
	{
		int lc = ggrid.lc[ iFace ];
		int rc = ggrid.rc[ iFace ];
		int lZone = cellzone[ lc ];
		int rZone = cellzone[ rc ];

		int ftype = ggrid.fTypes[ iFace ];

		if ( lZone != rZone )
		{
			int localCell_L = localCells[ lc ];
			int localCell_R = localCells[ rc ]; //Local cell count in another zone

			int bctype = -1;

			ScalarGrid & gridL = *grids[ lZone ];
			ScalarGrid & gridR = *grids[ rZone ];

			gridL.AddFaceType( ftype );
			gridR.AddFaceType( ftype );

			gridL.AddInterfaceBcFace( iFace, bctype, localCell_L, ONEFLOW::INVALID_INDEX, rZone, localCell_R );
			gridR.AddInterfaceBcFace( iFace, bctype, ONEFLOW::INVALID_INDEX, localCell_R, lZone, localCell_L );
		}
	}

	//Finally, scan the inner face of the block
	for ( int iFace = nBFaces; iFace < nFaces; ++ iFace )
	{
		int lc = ggrid.lc[ iFace ];
		int rc = ggrid.rc[ iFace ];
		int lZone = cellzone[ lc ];
		int rZone = cellzone[ rc ];

		if ( lZone == rZone )
		{
			//inner face bctype = 0
			int bctype = 0;
			int localCell_L = localCells[ lc ];
			int localCell_R = localCells[ rc ];

			ScalarGrid & grid = *grids[ lZone ];

			int ftype = ggrid.fTypes[ iFace ];

			grid.AddFaceType( ftype );
			grid.AddInnerFace( iFace, bctype, localCell_L, localCell_R );
		}
	}
    return grids;
}

void GridPartition::ReconstructInterfaceTopo( std::vector< std::unique_ptr< ScalarGrid > > & grids )
{
	int nZones = static_cast< int >( grids.size() );
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		if ( ! grids[ iZone ] || ! grids[ iZone ]->scalarIFace )
		{
			throw std::runtime_error( "GridPartition::ReconstructInterfaceTopo: zone has no interface topology" );
		}
		ScalarGrid & grid = *grids[ iZone ];
		ScalarIFace & scalarIFace = *grid.scalarIFace;
		int nNeis = static_cast< int >( scalarIFace.data.size() );
		for ( int iNei = 0; iNei < nNeis; ++ iNei )
		{
			ScalarIFaceIJ & iFaceIJ = scalarIFace.data[ iNei ];
			const int jZone = iFaceIJ.zonej;
			if ( jZone < 0 || jZone >= nZones || ! grids[ jZone ] || ! grids[ jZone ]->scalarIFace )
			{
				throw std::runtime_error( "GridPartition::ReconstructInterfaceTopo: interface references an invalid neighbor zone" );
			}
			if ( iFaceIJ.iglobalfaces.size() != iFaceIJ.ifaces.size() ||
				 iFaceIJ.iglobalfaces.size() != iFaceIJ.cells.size() )
			{
				throw std::runtime_error( "GridPartition::ReconstructInterfaceTopo: neighbor interface arrays have inconsistent sizes" );
			}
			std::cout << " iZone = " << iZone << " iNei = " << iNei << " jZone = " << jZone << "\n";
			grids[ jZone ]->scalarIFace->CalcLocalInterfaceId( iZone, iFaceIJ.iglobalfaces, iFaceIJ.target_ifaces );
		}
	}

	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		ScalarIFace & scalarIFace = *grids[ iZone ]->scalarIFace;
		const size_t nIFaces = scalarIFace.iglobalfaces.size();
		if ( scalarIFace.zones.size() != nIFaces || scalarIFace.cells.size() != nIFaces )
		{
			throw std::runtime_error( "GridPartition::ReconstructInterfaceTopo: interface mapping arrays have inconsistent sizes" );
		}
		for ( size_t iFace = 0; iFace < nIFaces; ++ iFace )
		{
			const int igface = scalarIFace.iglobalfaces[ iFace ];
			const int jZone = scalarIFace.zones[ iFace ];
			if ( jZone < 0 || jZone >= nZones || ! grids[ jZone ] || ! grids[ jZone ]->scalarIFace )
			{
				throw std::runtime_error( "GridPartition::ReconstructInterfaceTopo: interface mapping references an invalid neighbor zone" );
			}
			const int jlocalface = grids[ jZone ]->scalarIFace->GetLocalInterfaceId( igface );
			scalarIFace.target_interfaces.push_back( jlocalface );
		}
	}

}

void GridPartition::CalcInterfaceToBcFace( std::vector< std::unique_ptr< ScalarGrid > > & grids )
{
	int nZones = static_cast< int >( grids.size() );
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		grids[ iZone ]->CalcInterfaceToBcFace();
	}
}

void GridPartition::ReconstructNeighbor( std::vector< std::unique_ptr< ScalarGrid > > & grids )
{
	int nZones = static_cast< int >( grids.size() );
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		ScalarGrid & grid = *grids[ iZone ];
		ScalarIFace & scalarIFace = *grid.scalarIFace;
		scalarIFace.ReconstructNeighbor();
	}
}

void GridPartition::ReconstructNode( const ScalarGrid & ggrid, std::vector< std::unique_ptr< ScalarGrid > > & grids )
{
	int nZones = static_cast< int >( grids.size() );
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		ScalarGrid & grid = *grids[ iZone ];
		grid.ReconstructNode( ggrid );
		grid.Normalize();
		grid.CalcMetrics1D();
	}
}

EndNameSpace

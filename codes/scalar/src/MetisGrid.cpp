/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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
#include "HXSort.h"
#include "HXMath.h"
#include "Boundary.h"
#include <iostream>
#include <vector>
#include <algorithm>


BeginNameSpace( ONEFLOW )

MetisSplit::MetisSplit()
{
}

MetisSplit::~MetisSplit()
{
	;
}

void MetisSplit::ManualPartition( ScalarGrid * ggrid, int nPart, MetisIntList & cellzone )
{
	int nFaces = ggrid->GetNFaces();
	int nCells = ggrid->GetNCells();
	int nBFaces = ggrid->GetNBFaces();
	int nInnerFaces = nFaces - nBFaces;

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
}

void MetisSplit::MetisPartition( ScalarGrid * ggrid, int nPart, MetisIntList & cellzone )
{
	int nFaces = ggrid->GetNFaces();
	int nCells = ggrid->GetNCells();
	int nBFaces = ggrid->GetNBFaces();
	int nInnerFaces = nFaces - nBFaces;

	xadj.resize( nCells + 1 );
	adjncy.resize( 2 * nInnerFaces );
	cellzone.resize( nCells );

	if ( nPart == nCells )
	{
		ManualPartition( ggrid, nPart, cellzone );
		return;
	}

	ScalarGetXadjAdjncy( ggrid, xadj, adjncy );
	ScalarPartitionByMetis( nCells, xadj, adjncy, nPart, cellzone );
}

void MetisSplit::ScalarGetXadjAdjncy( ScalarGrid * ggrid, MetisIntList & xadj, MetisIntList & adjncy )
{   
	int nCells = ggrid->GetNCells();

	EList c2c;
	ggrid->CalcC2C( c2c );

	xadj[ 0 ]  = 0;
	int iCount = 0;
	for ( int iCell = 0; iCell < nCells; ++ iCell )
	{
		xadj[ iCell + 1 ] = xadj[ iCell ] + c2c[ iCell ].size();
		for ( int j = 0; j < c2c[ iCell ].size(); ++ j )
		{
			adjncy[ iCount ++ ] = c2c[ iCell ][ j ];
		}
	}
}

void MetisSplit::ScalarPartitionByMetis( idx_t nCells, MetisIntList & xadj, MetisIntList & adjncy, int nPart, MetisIntList & cellzone )
{
	idx_t   ncon     = 1;
	idx_t   * vwgt   = 0;
	idx_t   * vsize  = 0;
	idx_t   * adjwgt = 0;
	float * tpwgts = 0;
	float * ubvec  = 0;
	idx_t options[ METIS_NOPTIONS ];
	idx_t wgtflag = 0;
	idx_t numflag = 0;
	idx_t objval;
	idx_t nZone = nPart;

	METIS_SetDefaultOptions( options );
	std::cout << "Now begining partition graph!\n";
	if ( nZone > 8 )
	{
		std::cout << "Using K-way Partitioning!\n";
		METIS_PartGraphKway( & nCells, & ncon, & xadj[ 0 ], & adjncy[ 0 ], vwgt, vsize, adjwgt, 
			& nZone, tpwgts, ubvec, options, & objval, & cellzone[ 0 ] );
	}
	else
	{
		std::cout << "Using Recursive Partitioning!\n";
		METIS_PartGraphRecursive( & nCells, & ncon, & xadj[ 0 ], & adjncy[ 0 ], vwgt, vsize, adjwgt, 
			& nZone, tpwgts, ubvec, options, & objval, & cellzone[ 0 ] );
	}
	std::cout << "The interface number: " << objval << std::endl; 
	std::cout << "Partition is finished!\n";
}

GridPartition::GridPartition()
{
}

GridPartition::~GridPartition()
{
	;
}

void GridPartition::PartitionGrid( ScalarGrid * ggrid, int nPart, std::vector< ScalarGrid * > *grids )
{
	this->ggrid = ggrid;
	this->nPart = nPart;
	this->grids = grids;

	this->ReconstructGridFaceTopo();
	this->ReconstructNeighbor();
	this->ReconstructInterfaceTopo();
	this->CalcInterfaceToBcFace();
	this->ReconstructNode();
}

int GridPartition::GetNZones()
{
	return this->nPart;
}

void GridPartition::AllocateGrid( int nZones )
{
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		ScalarGrid * grid = new ScalarGrid();
		grid->id = iZone;
		( * this->grids ).push_back( grid );
	}
}

void GridPartition::ReconstructGridFaceTopo()
{
	//calc cellzone;
	MetisSplit metisSplit;
	MetisIntList cellzone;
	metisSplit.MetisPartition( this->ggrid, this->nPart, cellzone );

	this->AllocateGrid( this->nPart );

	int nZones = this->GetNZones();
	int nFaces = ggrid->GetNFaces();
	int nCells = ggrid->GetNCells();
	int nBFaces = ggrid->GetNBFaces();

	std::vector<int> zoneCount( nZones, 0 );
	std::vector<int> localCells; //global cell id -> local cell id
	localCells.resize( nCells );

	for ( int iCell = 0; iCell < nCells; ++ iCell )
	{
		int iZone = cellzone[ iCell ];
		localCells[ iCell ] = zoneCount[ iZone ] ++;
		int eType = this->ggrid->eTypes[ iCell ];
		ScalarGrid * grid = ( * this->grids )[ iZone ];
		grid->eTypes.AddData( eType );
	}

	//First scan the global physical boundary
	for ( int iFace = 0; iFace < nBFaces; ++ iFace )
	{
		int lc = ggrid->lc[ iFace ];
		int bctype = ggrid->bcTypes[ iFace ];
		int lZone = cellzone[ lc ];
		//global face id = iFace, local face id = faceid.size();
		//global face node: 20,10,30,40, local face node 1 2 4 3 for example
		//global coor x[20],y[20],z[20],x[10],y[10],z[10]
		//local coor x[1],y[1],z[1],x[2],y[2],z[2]
		int localCell = localCells[ lc ];
		int ftype = ggrid->fTypes[ iFace ];
		ScalarGrid * gridL = ( * this->grids )[ lZone ];
		gridL->AddFaceType( ftype );
		gridL->AddPhysicalBcFace( iFace, bctype, localCell, ONEFLOW::INVALID_INDEX );
	}

	//Then scan the internal block interface
	for ( int iFace = nBFaces; iFace < nFaces; ++ iFace )
	{
		int lc = ggrid->lc[ iFace ];
		int rc = ggrid->rc[ iFace ];
		int lZone = cellzone[ lc ];
		int rZone = cellzone[ rc ];

		int ftype = ggrid->fTypes[ iFace ];

		if ( lZone != rZone )
		{
			int localCell_L = localCells[ lc ];
			int localCell_R = localCells[ rc ]; //Local cell count in another zone

			int bctype = -1;

			ScalarGrid * gridL = ( * this->grids )[ lZone ];
			ScalarGrid * gridR = ( * this->grids )[ rZone ];

			gridL->AddFaceType( ftype );
			gridR->AddFaceType( ftype );

			gridL->AddInterfaceBcFace( iFace, bctype, localCell_L, ONEFLOW::INVALID_INDEX, rZone, localCell_R );
			gridR->AddInterfaceBcFace( iFace, bctype, ONEFLOW::INVALID_INDEX, localCell_R, lZone, localCell_L );
		}
	}

	//Finally, scan the inner face of the block
	for ( int iFace = nBFaces; iFace < nFaces; ++ iFace )
	{
		int lc = ggrid->lc[ iFace ];
		int rc = ggrid->rc[ iFace ];
		int lZone = cellzone[ lc ];
		int rZone = cellzone[ rc ];

		if ( lZone == rZone )
		{
			//inner face bctype = 0
			int bctype = 0;
			int localCell_L = localCells[ lc ];
			int localCell_R = localCells[ rc ];

			ScalarGrid * grid = ( * this->grids )[ lZone ];

			int ftype = ggrid->fTypes[ iFace ];

			grid->AddFaceType( ftype );
			grid->AddInnerFace( iFace, bctype, localCell_L, localCell_R );
		}
	}
}

void GridPartition::ReconstructInterfaceTopo()
{
	int nZones = this->GetNZones();
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		int nNeis = ( * this->grids )[ iZone ]->scalarIFace->data.size();
		for ( int iNei = 0; iNei < nNeis; ++ iNei )
		{
			ScalarIFaceIJ & iFaceIJ = ( * this->grids )[ iZone ]->scalarIFace->data[ iNei ];
			int jZone =  iFaceIJ.zonej;
			std::cout << " iZone = " << iZone << " iNei = " << iNei << " jZone = " << jZone << "\n";
			( * this->grids )[ jZone ]->scalarIFace->CalcLocalInterfaceId( iZone, iFaceIJ.iglobalfaces, iFaceIJ.target_ifaces );
		}
	}

	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		ScalarIFace * scalarIFace = ( * this->grids )[ iZone ]->scalarIFace;
		int nIFaces = scalarIFace->iglobalfaces.size();
		for ( int iFace = 0; iFace < nIFaces; ++ iFace )
		{
			int igface = scalarIFace->iglobalfaces[ iFace ];
			int jZone  = scalarIFace->zones[ iFace ];
			int jlocalface = ( * this->grids )[ jZone ]->scalarIFace->GetLocalInterfaceId( igface );
			scalarIFace->target_interfaces.push_back( jlocalface );
		}
	}

}

void GridPartition::CalcInterfaceToBcFace()
{
	int nZones = this->GetNZones();
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		( * this->grids )[ iZone ]->CalcInterfaceToBcFace();
	}
}

void GridPartition::ReconstructNeighbor()
{
	int nZones = this->GetNZones();
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		( * this->grids )[ iZone ]->scalarIFace->ReconstructNeighbor();
	}
}

void GridPartition::ReconstructNode()
{
	int nZones = this->GetNZones();
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		ScalarGrid * grid = ( * this->grids )[ iZone ];
		grid->ReconstructNode( ggrid );
		grid->Normalize();
		grid->CalcMetrics1D();
	}
}

EndNameSpace

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

#include "MetisGrid.h"
#include "ScalarGrid.h"
#include "Constant.h"
#include "HXCgns.h"
#include "ElementHome.h"
#include "HXSort.h"
#include "HXMath.h"
#include "Boundary.h"
#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

BeginNameSpace( ONEFLOW )

MetisPart::MetisPart( ScalarGrid * ggrid )
{
	this->ggrid = ggrid;
}

MetisPart::~MetisPart()
{
	;
}

void MetisPart::MetisPartition( int nPart, MetisIntList & cellzone )
{
	int nFaces = ggrid->GetNFaces();
	int nCells = ggrid->GetNCells();
	int nBFaces = ggrid->GetNBFaces();
	int nInnerFaces = nFaces - nBFaces;

	xadj.resize( nCells + 1 );
	adjncy.resize( 2 * nInnerFaces );
	cellzone.resize( nCells );

	ScalarGetXadjAdjncy( ggrid, xadj, adjncy );
	ScalarPartitionByMetis( nCells, xadj, adjncy, nPart, cellzone );
}

void MetisPart::ScalarGetXadjAdjncy( ScalarGrid * ggrid, MetisIntList & xadj, MetisIntList & adjncy )
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

void MetisPart::ScalarPartitionByMetis( idx_t nCells, MetisIntList & xadj, MetisIntList & adjncy, int nPart, MetisIntList & cellzone )
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
	cout << "Now begining partition graph!\n";
	if ( nZone > 8 )
	{
		cout << "Using K-way Partitioning!\n";
		METIS_PartGraphKway( & nCells, & ncon, & xadj[ 0 ], & adjncy[ 0 ], vwgt, vsize, adjwgt, 
			& nZone, tpwgts, ubvec, options, & objval, & cellzone[ 0 ] );
	}
	else
	{
		cout << "Using Recursive Partitioning!\n";
		METIS_PartGraphRecursive( & nCells, & ncon, & xadj[ 0 ], & adjncy[ 0 ], vwgt, vsize, adjwgt, 
			& nZone, tpwgts, ubvec, options, & objval, & cellzone[ 0 ] );
	}
	cout << "The interface number: " << objval << endl; 
	cout << "Partition is finished!\n";
}

NetGrid::NetGrid()
{
	;
}

NetGrid::~NetGrid()
{
	this->DeAllocateGrid();
}

void NetGrid::AllocateGrid( int nZones )
{
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		ScalarGrid * grid = new ScalarGrid();
		grids.push_back( grid );
	}
}

void NetGrid::DeAllocateGrid()
{
	int nZones = grids.size();
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		delete grids[ iZone ];
	}
}

GridTopo::GridTopo()
{
}

GridTopo::~GridTopo()
{
}

void GridTopo::AddFaceId( int iFace )
{
	faceid.push_back( iFace );
}

void GridTopo::AddFaceType( int faceType )
{
	this->facetype.push_back( faceType );
}

void GridTopo::ReconstructNode( EList & faces )
{
	int nFaces = faceid.size();
	for ( int iFace = 0; iFace < nFaces; ++ iFace )
	{
		int iGFace = faceid[ iFace ];
		vector< int > & face = faces[ iGFace ];
		int nNode = face.size();
		for ( int iNode = 0; iNode < nNode; ++ iNode )
		{
			nodeset.insert( face[ iNode ] );
		}
	}
}

GridTopos::GridTopos()
{
}

GridTopos::~GridTopos()
{
}

void GridTopos::Allocate( int nZones )
{
	data.resize( nZones );
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		data[ iZone ].zoneid = iZone;
	}
}

void GridTopos::CalcGlobal2LocalCells( MetisIntList & cellzone )
{
	int nZones = this->data.size();
	int nCells = cellzone.size();
	vector<int> zoneCount( nZones, 0 );
	gLCells.resize( nCells );

	for ( int iCell = 0; iCell < nCells; ++ iCell )
	{
		int iZone = cellzone[ iCell ];
		gLCells[ iCell ] = zoneCount[ iZone ] ++;
	}
}

void GridTopos::CalcInterface()
{
	int nZones = this->data.size();
	this->scalarIFaces.data.resize( nZones );
	this->scalarIFaces.CalcInterface( this );
	int kkk = 1;
}

Part::Part()
{
}

Part::~Part()
{
	;
}

void Part::PartitionGrid( ScalarGrid * ggrid, int nPart, NetGrid * netGrid )
{
	this->ggrid = ggrid;
	this->nPart = nPart;
	this->netGrid = netGrid;
	//calc cellzone;
	this->CalcCellZone();
	this->ReconstructAllZones();
}

void Part::CalcCellZone()
{
	MetisPart metisPart( this->ggrid );
	metisPart.MetisPartition( this->nPart, this->cellzone );
}

void Part::ReconstructAllZones()
{
	netGrid->AllocateGrid( this->nPart );
	this->CalcGlobalInterface();
	gtopos.Allocate( this->nPart );
	gtopos.CalcGlobal2LocalCells( this->cellzone );
	this->ReconstructGrid();
	this->ReconstructNode();
	gtopos.CalcInterface();
	int kkk = 1;
}

void Part::CalcGlobalInterface()
{
	int nFaces = ggrid->GetNFaces();
	int nCells = ggrid->GetNCells();
	int nBFaces = ggrid->GetNBFaces();
	int nInnerFaces = nFaces - nBFaces;
	//global interfaces
	vector< int > interfaces;
    for ( int iFace = nBFaces; iFace < nFaces; ++ iFace )
    {
		int lc = ggrid->lc[ iFace ];
		int rc = ggrid->rc[ iFace ];
		if ( cellzone[ lc ] != cellzone[ rc ] )
		{
			interfaces.push_back( iFace );
		}
    }
	int kkk = 1;
}

void Part::ReconstructGrid()
{
	int nFaces = ggrid->GetNFaces();
	int nCells = ggrid->GetNCells();
	int nBFaces = ggrid->GetNBFaces();

	for ( int iFace = 0; iFace < nBFaces; ++ iFace )
	{
		int lc = ggrid->lc[ iFace ];
		int lZone = cellzone[ lc ];
		gtopos[ lZone ].AddFaceId( iFace );
		gtopos[ lZone ].AddFaceType( 1 );
		int localCell = gtopos.gLCells[ lc ];
		gtopos[ lZone ].lc.AddData( localCell );
		gtopos[ lZone ].rc.AddData( ONEFLOW::INVALID_INDEX );
	}

	for ( int iFace = nBFaces; iFace < nFaces; ++ iFace )
	{
		int lc = ggrid->lc[ iFace ];
		int rc = ggrid->rc[ iFace ];
		int lZone = cellzone[ lc ];
		int rZone = cellzone[ rc ];

		if ( lZone == rZone )
		{
			gtopos[ lZone ].AddFaceId( iFace );
			gtopos[ lZone ].AddFaceType( 0 );

			int localCell_L = gtopos.gLCells[ lc ];
			int localCell_R = gtopos.gLCells[ rc ];

			gtopos[ lZone ].lc.AddData( localCell_L );
			gtopos[ lZone ].rc.AddData( localCell_R );
		}
		else
		{
			gtopos[ lZone ].AddFaceId( iFace );
			gtopos[ rZone ].AddFaceId( iFace );
			gtopos[ lZone ].AddFaceType( -1 );
			gtopos[ rZone ].AddFaceType( -1 );

			int localCell_L = gtopos.gLCells[ lc ];
			int localCell_R = gtopos.gLCells[ rc ]; //Local cell count in another zone

			gtopos[ lZone ].lc.AddData( localCell_L );
			gtopos[ lZone ].rc.AddData( ONEFLOW::INVALID_INDEX );

			gtopos[ rZone ].lc.AddData( ONEFLOW::INVALID_INDEX );
			gtopos[ rZone ].rc.AddData( localCell_R );
		}
	}
}

void Part::ReconstructNode()
{
	int nZones = this->nPart;
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		gtopos[ iZone ].ReconstructNode( ggrid->faces );
	}
}

EndNameSpace
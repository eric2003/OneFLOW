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

void MetisPart::ManualPartition( int nPart, MetisIntList & cellzone )
{
	int nFaces = ggrid->GetNFaces();
	int nCells = ggrid->GetNCells();
	int nBFaces = ggrid->GetNBFaces();
	int nInnerFaces = nFaces - nBFaces;

	vector< int > tmp;
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

void MetisPart::MetisPartition( int nPart, MetisIntList & cellzone )
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
		ManualPartition( nPart, cellzone );
		return;
	}

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

Part::Part()
{
}

Part::~Part()
{
	;
}

void Part::PartitionGrid( ScalarGrid * ggrid, int nPart, vector< ScalarGrid * > *grids )
{
	this->ggrid = ggrid;
	this->nPart = nPart;
	this->grids = grids;
	//calc cellzone;
	this->CalcCellZone();
	this->ReconstructAllZones();
}

void Part::CalcCellZone()
{
	MetisPart metisPart( this->ggrid );
	metisPart.MetisPartition( this->nPart, this->cellzone );
}

int Part::GetNZones()
{
	return this->nPart;
}

void Part::AllocateGrid( int nZones )
{
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		ScalarGrid * grid = new ScalarGrid( iZone );
		( * this->grids ).push_back( grid );
	}
}

void Part::ReconstructAllZones()
{
	this->AllocateGrid( this->nPart );
	this->ReconstructGridFaceTopo();
	this->ReconstructNeighbor();
	this->ReconstructInterfaceTopo();
	this->CalcInterfaceToBcFace();
	this->ReconstructNode();
}

void Part::ReconstructGridFaceTopo()
{
	int nZones = this->GetNZones();
	int nFaces = ggrid->GetNFaces();
	int nCells = ggrid->GetNCells();
	int nBFaces = ggrid->GetNBFaces();

	vector<int> zoneCount( nZones, 0 );
	gLCells.resize( nCells );

	for ( int iCell = 0; iCell < nCells; ++ iCell )
	{
		int iZone = cellzone[ iCell ];
		gLCells[ iCell ] = zoneCount[ iZone ] ++;
		int eType = this->ggrid->eTypes[ iCell ];
		( * this->grids )[ iZone ]->eTypes.AddData( eType );
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
		int localCell = this->gLCells[ lc ];
		( * this->grids )[ lZone ]->AddPhysicalBcFace( iFace, bctype, localCell, ONEFLOW::INVALID_INDEX );
	}

	//Then scan the internal block interface
	for ( int iFace = nBFaces; iFace < nFaces; ++ iFace )
	{
		int lc = ggrid->lc[ iFace ];
		int rc = ggrid->rc[ iFace ];
		int lZone = cellzone[ lc ];
		int rZone = cellzone[ rc ];

		if ( lZone != rZone )
		{
			int localCell_L = this->gLCells[ lc ];
			int localCell_R = this->gLCells[ rc ]; //Local cell count in another zone

			int bctype = -1;

			( * this->grids )[ lZone ]->AddInterfaceBcFace( iFace, bctype, localCell_L, ONEFLOW::INVALID_INDEX, rZone, localCell_R );
			( * this->grids )[ rZone ]->AddInterfaceBcFace( iFace, bctype, ONEFLOW::INVALID_INDEX, localCell_R, lZone, localCell_L );
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
			int localCell_L = this->gLCells[ lc ];
			int localCell_R = this->gLCells[ rc ];

			( * this->grids )[ lZone ]->AddInnerFace( iFace, bctype, localCell_L, localCell_R );
		}
	}
}

void Part::ReconstructInterfaceTopo()
{
	int nZones = this->GetNZones();
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		int nNeis = ( * this->grids )[ iZone ]->scalarIFace->data.size();
		for ( int iNei = 0; iNei < nNeis; ++ iNei )
		{
			ScalarIFaceIJ & iFaceIJ = ( * this->grids )[ iZone ]->scalarIFace->data[ iNei ];
			int jZone =  iFaceIJ.zonej;
			//iZone的第iNei个邻居为jZone,iZone和jZone的交界面的global interface id为：
			//iFaceIJ.iglobalfaces,iFaceIJ.target_ifaces是这些interface在jZone里面的局部id
			//这些id由jZone计算发送给iZone的里iNei个信息存储
			//实际上这个信息iZone用不到，是jZone接收时使用的。
			cout << " iZone = " << iZone << " iNei = " << iNei << " jZone = " << jZone << "\n";
			( * this->grids )[ jZone ]->scalarIFace->CalcLocalInterfaceId( iZone, iFaceIJ.iglobalfaces, iFaceIJ.target_ifaces );
		}
	}
}

void Part::CalcInterfaceToBcFace()
{
	int nZones = this->GetNZones();
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		( * this->grids )[ iZone ]->CalcInterfaceToBcFace();
	}
}

void Part::ReconstructNeighbor()
{
	int nZones = this->GetNZones();
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		( * this->grids )[ iZone ]->scalarIFace->ReconstructNeighbor();
	}
}

void Part::ReconstructNode()
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
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

GridTopo::GridTopo( ScalarGrid * grid )
{
	this->grid = grid;
	this->scalarIFace = new ScalarIFace( grid->grid_id );
}

GridTopo::~GridTopo()
{
	delete this->scalarIFace;
}

int GridTopo::GetNBFaces()
{
	return this->bctypes.GetNElements();
}

void GridTopo::SetNCells( int nCells )
{
	this->nCells = nCells;
}

int GridTopo::GetNCells()
{
	return this->nCells;
}

void GridTopo::AddPhysicalBcFace( int global_face_id, int bctype, int lcell, int rcell )
{
	this->faceid.push_back( global_face_id );
	this->bctypes.AddData( bctype );
	this->facetype.AddData( bctype );

	this->grid->lc.AddData( lcell );
	this->grid->rc.AddData( rcell );
}

void GridTopo::AddInterfaceBcFace( int global_face_id, int bctype, int lcell, int rcell, int nei_zoneid, int nei_cellid )
{
	this->faceid.push_back( global_face_id );
	this->bctypes.AddData( bctype );
	this->facetype.AddData( bctype );

	this->grid->lc.AddData( lcell );
	this->grid->rc.AddData( rcell );

	this->AddInterface( global_face_id, nei_zoneid, nei_cellid );
}

void GridTopo::AddInnerFace( int global_face_id, int bctype, int lcell, int rcell )
{
	this->faceid.push_back( global_face_id );
	this->facetype.AddData( bctype );

	this->grid->lc.AddData( lcell );
	this->grid->rc.AddData( rcell );
}

void GridTopo::AddInterface( int global_interface_id, int neighbor_zoneid, int neighbor_cellid )
{
	this->scalarIFace->AddInterface( global_interface_id, neighbor_zoneid, neighbor_cellid );
}

void GridTopo::CalcInterfaceToBcFace()
{
	if ( scalarIFace->GetNIFaces() == 0 ) return;

	int nBFaces = this->GetNBFaces();

	scalarIFace->interface_to_bcface.resize( 0 );

	for ( int iBFace = 0; iBFace < nBFaces; ++ iBFace )
	{
		if ( ! BC::IsInterfaceBc( this->bctypes[ iBFace ] ) )
		{
			continue;
		}

		scalarIFace->interface_to_bcface.push_back( iBFace );
	}
}

void GridTopo::DumpGridInfo()
{
	cout << " DumpGridInfo Zone ID = " << this->grid->grid_id << "\n";
	int nIFace = this->scalarIFace->iglobalfaces.size();
	this->scalarIFace->DumpInterfaceMap();
}

void GridTopo::ReconstructNode( ScalarGrid * ggrid )
{
	this->ReconstructNode( ggrid->faces );
	this->CalcCoor( ggrid );
}

void GridTopo::CopyGrid( ScalarGrid * grid )
{
	grid->xn = this->xn;
	grid->yn = this->yn;
	grid->zn = this->zn;

	grid->faces = this->local_faces;

	grid->fBcTypes = this->bctypes;
	grid->eTypes.Resize( this->GetNCells() );

	grid->bcTypes = this->bctypes;

	this->Normalize( grid );

	grid->CalcMetrics1D();
}

void GridTopo::Normalize( ScalarGrid * grid )
{
	grid->faces = this->local_faces;

	int nFaces = grid->faces.GetNElements();
	for ( int iFace = 0; iFace < nFaces; ++ iFace )
	{
		if ( grid->lc[ iFace ] < 0 )
		{
			//need to reverse the node ordering
			vector< int > & face = grid->faces[ iFace ];
			std::reverse( face.begin(), face.end() );
			// now reverse lc and rc
			ONEFLOW::SWAP( grid->lc[ iFace ], grid->rc[ iFace ] );
		}
	}

	grid->SetBcGhostCell();
}

void GridTopo::ReconstructNode( EList & global_faces )
{
	int nFaces = faceid.size();
	for ( int iFace = 0; iFace < nFaces; ++ iFace )
	{
		int iGFace = faceid[ iFace ];
		vector< int > & face = global_faces[ iGFace ];
		int nNode = face.size();
		for ( int iNode = 0; iNode < nNode; ++ iNode )
		{
			nodeset.insert( face[ iNode ] );
		}
		local_faces.AddElem( face );
		faces.AddElem( face );
	}

	this->CalcGlobal2LocalNodeMap();
	this->CalcLocalFaceNodes();
}

void GridTopo::CalcGlobal2LocalNodeMap()
{
	int count = 0;
	for ( set<int>::iterator iter = nodeset.begin(); iter != nodeset.end(); ++ iter )
	{
		global_local_node.insert( pair<int, int>( *iter, count ++ ) );
	}
	int kkk = 1;
}

void GridTopo::CalcLocalFaceNodes()
{
	//local_faces
	int nFaces = local_faces.GetNElements();
	for ( int iFace = 0; iFace < nFaces; ++ iFace )
	{
		vector< int > & face = local_faces[ iFace ];
		int nNode = face.size();
		for ( int iNode = 0; iNode < nNode; ++ iNode )
		{
			int glbal_node_id = face[ iNode ];
			face[ iNode ] = global_local_node[ glbal_node_id ];
		}
	}
	int kkk = 1;
}

void GridTopo::ReorderInterface()
{
	int nFaces = faceid.size();
	for ( int iFace = 0; iFace < nFaces; ++ iFace )
	{
	}
}

void GridTopo::ReconstructNeighbor()
{
	this->scalarIFace->ReconstructNeighbor();
}

void GridTopo::CalcCoor( ScalarGrid * grid )
{
	for ( set<int>::iterator iter = nodeset.begin(); iter != nodeset.end(); ++ iter )
	{
		int iNode = *iter;
		Real xm = grid->xn[ iNode ];
		Real ym = grid->yn[ iNode ];
		Real zm = grid->zn[ iNode ];

		this->xn.AddData( xm );
		this->yn.AddData( ym );
		this->zn.AddData( zm );
	}
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

void Part::CalcGlobal2LocalCells( MetisIntList & cellzone )
{
	int nZones = this->GetNZones();
	int nCells = cellzone.size();
	vector<int> zoneCount( nZones, 0 );
	gLCells.resize( nCells );

	for ( int iCell = 0; iCell < nCells; ++ iCell )
	{
		int iZone = cellzone[ iCell ];
		gLCells[ iCell ] = zoneCount[ iZone ] ++;
	}

	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		( * this->grids )[ iZone ]->gridTopo->SetNCells( zoneCount[ iZone ] );
	}
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
	this->CalcGlobal2LocalCells( this->cellzone );
	this->ReconstructGridFaceTopo();
	this->ReconstructNeighbor();
	this->ReconstructInterfaceTopo();
	this->CalcInterfaceToBcFace();
	//this->DumpGridInfo();
	this->ReconstructNode();
	int kkk = 1;
}

void Part::DumpGridInfo()
{
	int nZones = this->GetNZones();
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		( * this->grids )[ iZone ]->gridTopo->DumpGridInfo();
	}
}

void Part::ReconstructGridFaceTopo()
{
	int nFaces = ggrid->GetNFaces();
	int nCells = ggrid->GetNCells();
	int nBFaces = ggrid->GetNBFaces();

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
		( * this->grids )[ lZone ]->gridTopo->AddPhysicalBcFace( iFace, bctype, localCell, ONEFLOW::INVALID_INDEX );
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

			( * this->grids )[ lZone ]->gridTopo->AddInterfaceBcFace( iFace, bctype, localCell_L, ONEFLOW::INVALID_INDEX, rZone, localCell_R );
			( * this->grids )[ rZone ]->gridTopo->AddInterfaceBcFace( iFace, bctype, ONEFLOW::INVALID_INDEX, localCell_R, lZone, localCell_L );
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

			( * this->grids )[ lZone ]->gridTopo->AddInnerFace( iFace, bctype, localCell_L, localCell_R );
		}
	}
}

void Part::ReconstructInterfaceTopo()
{
	int nZones = this->GetNZones();
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		int nNeis = ( * this->grids )[ iZone ]->gridTopo->scalarIFace->data.size();
		for ( int iNei = 0; iNei < nNeis; ++ iNei )
		{
			ScalarIFaceIJ & iFaceIJ = ( * this->grids )[ iZone ]->gridTopo->scalarIFace->data[ iNei ];
			int jZone =  iFaceIJ.zonej;
			//iZone的第iNei个邻居为jZone,iZone和jZone的交界面的global interface id为：
			//iFaceIJ.iglobalfaces,iFaceIJ.target_ifaces是这些interface在jZone里面的局部id
			//这些id由jZone计算发送给iZone的里iNei个信息存储
			//实际上这个信息iZone用不到，是jZone接收时使用的。
			cout << " iZone = " << iZone << " iNei = " << iNei << " jZone = " << jZone << "\n";
			( * this->grids )[ jZone ]->gridTopo->scalarIFace->CalcLocalInterfaceId( iZone, iFaceIJ.iglobalfaces, iFaceIJ.target_ifaces );
		}
	}
}

void Part::CalcInterfaceToBcFace()
{
	int nZones = this->GetNZones();
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		( * this->grids )[ iZone ]->gridTopo->CalcInterfaceToBcFace();
	}
}

void Part::ReconstructNeighbor()
{
	int nZones = this->GetNZones();
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		( * this->grids )[ iZone ]->gridTopo->ReconstructNeighbor();
	}
}

void Part::ReconstructNode()
{
	int nZones = this->GetNZones();
	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		( * this->grids )[ iZone ]->gridTopo->ReconstructNode( ggrid );
	}

	for ( int iZone = 0; iZone < nZones; ++ iZone )
	{
		( * this->grids )[ iZone ]->gridTopo->CopyGrid( ( * this->grids )[ iZone ] );
	}
}

EndNameSpace
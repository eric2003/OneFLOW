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

#include "ScalarGrid.h"
#include "DataBase.h"
#include <memory>
#include "ScalarCgns.h"
#include "CgnsZsection.h"
#include "CgnsSection.h"
#include "CgnsZbc.h"
#include "CgnsZbcBoco.h"
#include "CgnsBcBoco.h"
#include "CgnsCoor.h"
#include "HXStd.h"
#include "NodeMesh.h"
#include "GridState.h"
#include "CgnsZbase.h"
#include "CgnsBase.h"
#include "CgnsZone.h"
#include "CgnsFile.h"
#include "Constant.h"
#include "HXCgns.h"
#include "StringUtils.h"
#include "Dimension.h"
#include "ElementHome.h"
#include "HXMath.h"
#include "HXMathExt.h"
#include "HXLookup.h"
#include "Boundary.h"
#include "MetisGrid.h"
#include "ScalarIFace.h"
#include "DataBook.h"
#include "DataBaseIO.h"
#include "Prj.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <utility>
#include <stdexcept>
#include <limits>
#include <unordered_set>


BeginNameSpace( ONEFLOW )

namespace
{
void ValidateReorderMap( const IntList & orderMap, size_t elementCount, const char * operation )
{
	if ( orderMap.GetNElements() != elementCount )
	{
		throw std::invalid_argument( std::string( operation ) +
			": order map size must match list size" );
	}

	std::vector< bool > seen( elementCount, false );
	for ( size_t i = 0; i < elementCount; ++ i )
	{
		const int sourceIndex = orderMap[ i ];
		if ( sourceIndex < 0 || static_cast< size_t >( sourceIndex ) >= elementCount ||
			 seen[ static_cast< size_t >( sourceIndex ) ] )
		{
			throw std::invalid_argument( std::string( operation ) +
				": order map must be a permutation" );
		}
		seen[ static_cast< size_t >( sourceIndex ) ] = true;
	}
}

template< typename T >
void ReserveForAppend( std::vector< T > & values )
{
	if ( values.size() < values.capacity() )
	{
		return;
	}

	const size_t requiredCapacity = values.size() + 1;
	const size_t maxCapacity = values.max_size();
	if ( requiredCapacity > maxCapacity )
	{
		throw std::length_error( "ScalarGrid: append exceeds vector maximum size" );
	}

	const size_t currentCapacity = values.capacity();
	size_t newCapacity = currentCapacity == 0 ? 1 : currentCapacity;
	if ( newCapacity <= maxCapacity / 2 )
	{
		newCapacity *= 2;
	}
	else
	{
		newCapacity = maxCapacity;
	}
	if ( newCapacity < requiredCapacity )
	{
		newCapacity = requiredCapacity;
	}
	values.reserve( newCapacity );
}
}

size_t RealList::GetNElements() const
{
	return data.size();
}

void RealList::AddData( Real value )
{
	data.push_back( value );
}

void RealList::Reserve( int new_capacity )
{
	if ( new_capacity < 0 )
	{
		throw std::invalid_argument( "RealList::Reserve: capacity must be non-negative" );
	}
	data.reserve( static_cast< size_t >( new_capacity ) );
}

void RealList::Resize( int new_size )
{
	if ( new_size < 0 )
	{
		throw std::invalid_argument( "RealList::Resize: size must be non-negative" );
	}
	data.resize( static_cast< size_t >( new_size ) );
}

IntList::IntList( const IntList & rhs )
{
	this->data = rhs.data;
}

size_t IntList::GetNElements() const
{
	return data.size();
}

void IntList::AddData( int value )
{
	data.push_back( value );
}

void IntList::Resize( int new_size )
{
	if ( new_size < 0 )
	{
		throw std::invalid_argument( "IntList::Resize: size must be non-negative" );
	}
	data.resize( static_cast< size_t >( new_size ) );
}

void IntList::Reserve( int new_size )
{
	if ( new_size < 0 )
	{
		throw std::invalid_argument( "IntList::Reserve: size must be non-negative" );
	}
	data.reserve( static_cast< size_t >( new_size ) );
}

void IntList::ReOrder( const IntList & orderMap )
{
	const size_t nElems = this->GetNElements();
	ValidateReorderMap( orderMap, nElems, "IntList::ReOrder" );

	const std::vector< int > original = data;
	for ( size_t i = 0; i < nElems; ++ i )
	{
		data[ i ] = original[ static_cast< size_t >( orderMap[ i ] ) ];
	}
}

size_t EList::GetNElements() const
{
	return data.size();
}

void EList::AddElem( const IntList &elem )
{
	this->data.push_back( elem.data );
}

void EList::AddElem( const std::vector< int > &elem )
{
	this->data.push_back( elem );
}

void EList::ReOrder( const IntList & orderMap )
{
	const size_t nElems = this->GetNElements();
	ValidateReorderMap( orderMap, nElems, "EList::ReOrder" );

	const std::vector< std::vector< int > > original = data;
	for ( size_t i = 0; i < nElems; ++ i )
	{
		data[ i ] = original[ static_cast< size_t >( orderMap[ i ] ) ];
	}
}

void EList::Resize( int new_size )
{
	if ( new_size < 0 )
	{
		throw std::invalid_argument( "EList::Resize: size must be non-negative" );
	}
	data.resize( static_cast< size_t >( new_size ) );
}

void EList::Reserve( int new_size )
{
	if ( new_size < 0 )
	{
		throw std::invalid_argument( "EList::Reserve: size must be non-negative" );
	}
	data.reserve( static_cast< size_t >( new_size ) );
}

void ScalarBcco::PushBoundaryFace( int pt, int eType )
{
	IntList elem;
	elem.AddData( pt );

	this->elements.AddElem( elem );
	this->eTypes.AddData( eType );
}

void ScalarBcco::AddBcPoint( int bcVertex )
{
	vertexList.push_back( bcVertex );
}

void ScalarBcco::ScanBcFace( ScalarGrid & grid )
{
	std::cout << " BCTypeName = " << ONEFLOW::GetCgnsBcName( this->bcType ) << std::endl;

	IntSet bcVertex;
	this->ProcessVertexBc( bcVertex );
	
	grid.ScanBcFace( bcVertex, this->bcType );
}

void ScalarBcco::ProcessVertexBc( IntSet & bcVertex )
{
	for ( int iBcPoint = 0; iBcPoint < this->vertexList.size(); ++ iBcPoint )
	{
		bcVertex.insert( this->vertexList[ iBcPoint ] );
	}
}


ScalarBccos::~ScalarBccos() = default;

void ScalarBccos::AddBcco( std::unique_ptr< ScalarBcco > scalarBcco )
{
	this->bccos.push_back( std::move( scalarBcco ) );
}

void ScalarBccos::ScanBcFace( ScalarGrid & grid )
{
	for ( int iBoco = 0; iBoco < this->bccos.size(); ++ iBoco )
	{
		std::cout << " iBoco = " << iBoco << " ";
		ScalarBcco & scalarBcco = *this->bccos[ iBoco ];
		scalarBcco.ScanBcFace( grid );
	}
}

ScalarGrid::ScalarGrid()
	: dataBase( std::make_unique< DataBase >() ),
	  scalarBccos( std::make_unique< ScalarBccos >() ),
	  scalarIFace( std::make_unique< ScalarIFace >() )
{
	this->nNodes = 0;
	this->nCells = 0;
	this->nBFaces = 0;
	this->nFaces = 0;
	this->nTCells = 0;
	this->id = 0;
	this->localId = 0;
	this->level = 0;
	this->volBcType = -1;
	this->type = ONEFLOW::UMESH;
}

ScalarGrid::~ScalarGrid() = default;

DataBase * ScalarGrid::GetDataBase()
{
	return dataBase.get();
}

const DataBase * ScalarGrid::GetDataBase() const
{
	return dataBase.get();
}

DataBase & ScalarGrid::RequireDataBase()
{
	if ( dataBase == nullptr )
	{
		throw std::logic_error( "ScalarGrid: DataBase is not initialized" );
	}
	return *dataBase;
}

const DataBase & ScalarGrid::RequireDataBase() const
{
	if ( dataBase == nullptr )
	{
		throw std::logic_error( "ScalarGrid: DataBase is not initialized" );
	}
	return *dataBase;
}

void ScalarGrid::ResetMeshData()
{
	xn.data.clear();
	yn.data.clear();
	zn.data.clear();

	// Reuse the topology reset so mesh and topology rebuilds share one lifecycle path.
	this->ResetTopologyData();

	elements.data.clear();
	boundaryElements.data.clear();

	bcETypes.data.clear();
	eTypes.data.clear();
	bcNameIds.data.clear();

	scalarBccos = std::make_unique< ScalarBccos >();
	scalarIFace = std::make_unique< ScalarIFace >();

	nNodes = 0;
	nCells = 0;
}

void ScalarGrid::ResetTopologyData()
{
	lc.data.clear();
	rc.data.clear();
	lpos.data.clear();
	rpos.data.clear();
	faces.data.clear();
	fTypes.data.clear();
	fBcTypes.data.clear();
	bcTypes.data.clear();

	// These mappings are derived from the face topology and cannot survive its reset.
	cell2faces.clear();
	c2fpos.clear();
	global_faceid.clear();

	// Interface IDs, neighbor groups, and boundary-face links all refer to this topology.
	if ( scalarIFace != nullptr )
	{
		scalarIFace->ResetTopology();
	}

	nFaces = 0;
	nBFaces = 0;

	// Geometry is derived from topology and must be invalidated with it.
	this->ResetGeometryData();
}

void ScalarGrid::ValidateFaceTopology( const char * operation ) const
{
	if ( this->nFaces < 0 || this->nBFaces < 0 || this->nBFaces > this->nFaces ||
		 this->nBFaces != this->GetNBFaces() )
	{
		throw std::logic_error( std::string( operation ) +
			": face and boundary-face counts are inconsistent" );
	}

	const size_t faceCount = static_cast< size_t >( this->nFaces );
	if ( this->faces.GetNElements() != faceCount ||
		 this->fTypes.GetNElements() != faceCount ||
		 this->lc.GetNElements() != faceCount ||
		 this->rc.GetNElements() != faceCount )
	{
		throw std::logic_error( std::string( operation ) +
			": face connectivity, types, and adjacency arrays must match nFaces" );
	}

	const int nodeCount = this->GetNNodes();
	const int cellCount = this->GetNCells();
	if ( this->nNodes != nodeCount || this->yn.GetNElements() != static_cast< size_t >( nodeCount ) ||
		 this->zn.GetNElements() != static_cast< size_t >( nodeCount ) )
	{
		throw std::logic_error( std::string( operation ) +
			": node count and coordinate arrays have inconsistent sizes" );
	}
	if ( this->nCells != cellCount )
	{
		throw std::logic_error( std::string( operation ) +
			": cell count and element type array have inconsistent sizes" );
	}

	for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
	{
		const std::vector< int > & face = this->faces[ iFace ];
		if ( face.empty() )
		{
			throw std::logic_error( std::string( operation ) +
				": face connectivity must not be empty" );
		}

		std::unordered_set< int > uniqueNodeIds;
		for ( const int nodeId : face )
		{
			if ( nodeId < 0 || nodeId >= nodeCount )
			{
				throw std::out_of_range( std::string( operation ) +
					": face references an invalid node index" );
			}
			if ( ! uniqueNodeIds.insert( nodeId ).second )
			{
				throw std::logic_error( std::string( operation ) +
					": face connectivity contains duplicate node indices" );
			}
		}

		const int leftCell = this->lc[ iFace ];
		const int rightCell = this->rc[ iFace ];
		if ( iFace < this->nBFaces )
		{
			const int expectedGhostCell = cellCount + iFace;
			const bool leftIsPhysical = leftCell >= 0 && leftCell < cellCount;
			const bool rightIsBoundarySentinel = rightCell == -1;
			const bool rightIsExpectedGhost = rightCell == expectedGhostCell;

			if ( ! leftIsPhysical ||
				 ( ! rightIsBoundarySentinel && ! rightIsExpectedGhost ) )
			{
				throw std::out_of_range( std::string( operation ) +
					": boundary face must have exactly one physical adjacent cell" );
			}
		}
		else
		{
			if ( leftCell < 0 || leftCell >= cellCount ||
				 rightCell < 0 || rightCell >= cellCount )
			{
				throw std::out_of_range( std::string( operation ) +
					": internal face references an invalid physical cell" );
			}
			if ( leftCell == rightCell )
			{
				throw std::logic_error( std::string( operation ) +
					": internal face must connect distinct cells" );
			}
		}
	}
}

void ScalarGrid::ResetGeometryData()
{
	xfc.data.clear();
	yfc.data.clear();
	zfc.data.clear();
	xfn.data.clear();
	yfn.data.clear();
	zfn.data.clear();
	area.data.clear();
	xcc.data.clear();
	ycc.data.clear();
	zcc.data.clear();
	vol.data.clear();

	nTCells = 0;
}

int ScalarGrid::GetNNodes() const
{
	return this->xn.GetNElements();
}

int ScalarGrid::GetNCells() const
{
	return this->eTypes.GetNElements();
}

int ScalarGrid::GetNTCells() const
{
	return this->GetNBFaces() + this->GetNCells();
}

int ScalarGrid::GetNFaces() const
{
	return this->faces.GetNElements();
}

int ScalarGrid::GetNBFaces() const
{
	return this->bcTypes.GetNElements();
}

void ScalarGrid::GenerateGrid( int ni, Real xmin, Real xmax )
{
	if ( ni < 2 )
	{
		throw std::invalid_argument( "ScalarGrid::GenerateGrid: ni must be at least 2" );
	}
	const Real domainLength = xmax - xmin;
	if ( ! std::isfinite( xmin ) || ! std::isfinite( xmax ) ||
		 ! std::isfinite( domainLength ) || domainLength <= 0.0 )
	{
		throw std::invalid_argument(
			"ScalarGrid::GenerateGrid: domain bounds must define a finite, positive interval" );
	}

	Real dx = domainLength / ( ni - 1 );
	if ( ! std::isfinite( dx ) || dx <= 0.0 )
	{
		throw std::invalid_argument(
			"ScalarGrid::GenerateGrid: spacing must be finite and positive" );
	}

	// A positive interval can still be too small to represent ni distinct
	// floating-point coordinates. Check representability before replacing the
	// current mesh, and keep the requested right endpoint exact.
	Real previousX = xmin;
	for ( int i = 1; i < ni; ++ i )
	{
		const Real x = i == ni - 1 ? xmax : xmin + i * dx;
		if ( ! std::isfinite( x ) || x <= previousX )
		{
			throw std::invalid_argument(
				"ScalarGrid::GenerateGrid: coordinates must be finite and strictly increasing" );
		}
		previousX = x;
	}

	// Validate scalar inputs and coordinate representability before clearing an existing mesh.
	this->ResetMeshData();

	// Reserve source-mesh storage once when the final sizes are known.
	// Topology is still built separately by CalcTopology().
	xn.Reserve( ni );
	yn.Reserve( ni );
	zn.Reserve( ni );
	elements.Reserve( ni - 1 );
	eTypes.Reserve( ni - 1 );

	for ( int i = 0; i < ni; ++ i )
	{
		Real xm = i == ni - 1 ? xmax : xmin + i * dx;
		Real ym = 0.0;
		Real zm = 0.0;

		xn.AddData( xm );
		yn.AddData( ym );
		zn.AddData( zm );
	}

	int ptL = 0;
	int ptR = ni - 1;

	for ( int i = 0; i < ni - 1; ++ i )
	{
		int p1 = i;
		int p2 = i + 1;

		const int eType = ONEFLOW::BAR_2;

		// This batch builder invalidated derived data once in ResetMeshData().
		// Append source connectivity directly instead of resetting topology per cell.
		this->elements.AddElem( std::vector< int >{ p1, p2 } );
		this->eTypes.AddData( eType );
	}

	auto scalarBccoL = std::make_unique< ScalarBcco >();
	scalarBccoL->bcName = "LeftOutFlow";
	scalarBccoL->bcType = ONEFLOW::BCOutflow;
	scalarBccoL->PushBoundaryFace( ptL, ONEFLOW::NODE );
	scalarBccoL->AddBcPoint( ptL );
	scalarBccos->AddBcco( std::move( scalarBccoL ) );

	auto scalarBccoR = std::make_unique< ScalarBcco >();
	scalarBccoR->bcName = "RightOutFlow";
	scalarBccoR->bcType = ONEFLOW::BCOutflow;
	scalarBccoR->PushBoundaryFace( ptR, ONEFLOW::NODE );
	scalarBccoR->AddBcPoint( ptR );
	scalarBccos->AddBcco( std::move( scalarBccoR ) );

}

void ScalarGrid::CalcVolumeSection( SectionManager & volumeSectionManager )
{
	IntSet typeSet;
	int nElements = this->elements.GetNElements();
	for ( int iElement = 0; iElement < nElements; ++ iElement )
	{
		int eType = this->eTypes[ iElement ];
		typeSet.insert( eType );
	}

	IntField cgns_types;

	ONEFLOW::Set2Array( typeSet, cgns_types );

	int nElementTypes = cgns_types.size();
	volumeSectionManager.Alloc( nElementTypes );

	for ( int iType = 0; iType < nElementTypes; ++ iType )
	{
		int current_eType = cgns_types[ iType ];
		SectionMarker & sectionMarker = *volumeSectionManager.data[ iType ];
		sectionMarker.cgns_type = current_eType;
		sectionMarker.name = ElementTypeName[ sectionMarker.cgns_type ];
		for ( int iElement = 0; iElement < nElements; ++ iElement )
		{
			int eType = this->eTypes[ iElement ];
			if ( eType == current_eType )
			{
				sectionMarker.elements.push_back( this->elements[ iElement ] );
				sectionMarker.elementIds.push_back( iElement );
			}
		}
		sectionMarker.nElements = sectionMarker.elements.size();
	}
}

void ScalarGrid::CalcBoundarySection( SectionManager & bcSectionManager )
{
	IntSet typeSet;

	int nBccos = scalarBccos->bccos.size();
	for ( int iBcco = 0; iBcco < nBccos; ++ iBcco )
	{
		ScalarBcco & scalarBcco = *scalarBccos->bccos[ iBcco ];
		int nElements = scalarBcco.eTypes.GetNElements();

		for ( int iElement = 0; iElement < nElements; ++ iElement )
		{
			int eType = scalarBcco.eTypes[ iElement ];
			typeSet.insert( eType );
		}
	}

	IntField cgns_types;

	ONEFLOW::Set2Array( typeSet, cgns_types );

	for ( int iBcco = 0; iBcco < nBccos; ++ iBcco )
	{
		ScalarBcco & scalarBcco = *scalarBccos->bccos[ iBcco ];
		int nElements = scalarBcco.eTypes.GetNElements();
		scalarBcco.local_globalIds.Resize( nElements );
	}

	int nElementTypes = cgns_types.size();
	bcSectionManager.Alloc( nElementTypes );
	int globalElementId = 0;
	for ( int iType = 0; iType < nElementTypes; ++ iType )
	{
		int current_eType = cgns_types[ iType ];
		SectionMarker & sectionMarker = *bcSectionManager.data[ iType ];
		sectionMarker.cgns_type = current_eType;
		sectionMarker.name = ElementTypeName[ sectionMarker.cgns_type ];
		for ( int iBcco = 0; iBcco < nBccos; ++ iBcco )
		{
			ScalarBcco & scalarBcco = *scalarBccos->bccos[ iBcco ];
			int nElements = scalarBcco.eTypes.GetNElements();

			for ( int iElement = 0; iElement < nElements; ++ iElement )
			{
				int eType = scalarBcco.eTypes[ iElement ];
				if ( eType == current_eType )
				{
					sectionMarker.elements.push_back( scalarBcco.elements[ iElement ] );
					sectionMarker.elementIds.push_back( globalElementId );
					scalarBcco.local_globalIds[ iElement ] = globalElementId;
					++ globalElementId;
				}
			}

		}

		sectionMarker.nElements = sectionMarker.elements.size();
	}
}

void ScalarGrid::SetCgnsZone( CgnsZone & cgnsZone )
{
	cgnsZone.zoneName = ONEFLOW::AddString( "Zone", cgnsZone.zId );
	cgnsZone.cgnsZoneType = CGNS_ENUMV( Unstructured );

	int nNodes = this->GetNNodes();
	int nCells = this->GetNCells();

	/* vertex size */
	cgnsZone.isize[ 0 ] = nNodes;
	/* cell size */
	cgnsZone.isize[ 1 ] = nCells;
	/* boundary vertex size (zero if elements not sorted) */
	cgnsZone.isize[ 2 ] = 0;

	CgnsCoor & cgnsCoor = cgnsZone.RequireCgnsCoor();

	cgnsCoor.SetNNode( nNodes );
	cgnsCoor.SetNCell( nCells );
	cgnsCoor.nCoor = 3;

	cgnsCoor.coorNameList[ 0 ] = "X";
	cgnsCoor.coorNameList[ 1 ] = "Y";
	cgnsCoor.coorNameList[ 2 ] = "Z";

	NodeMesh * nodeMesh = cgnsCoor.GetNodeMesh();
	nodeMesh->SetCoordinates( this->xn.data, this->yn.data, this->zn.data );

	DataType_t dataType = RealDouble;

	for ( int iCoor = 0; iCoor < cgnsCoor.nCoor; ++ iCoor )
	{
		int coordId = iCoor + 1;
		cgnsCoor.typeList[ iCoor ] = dataType;
		cgnsCoor.nNodeList[ iCoor ] = nNodes;
		cgnsCoor.Alloc( iCoor, static_cast<int>( nNodes ), dataType );
	}

	cgnsCoor.SetAllCoorData();

	SectionManager volSec;
	SectionManager bcSec;

	this->CalcVolumeSection( volSec );
	this->CalcBoundarySection( bcSec );

	int nVolSections = volSec.GetNSections();
	int nBcSections = bcSec.GetNSections();

	int nTotalSections = nVolSections + nBcSections;

	CgnsZsection & cgnsZsection = cgnsZone.RequireCgnsZsection();

	cgnsZsection.CreateCgnsSections( nTotalSections );

	int nVolCell = volSec.CalcTotalElem();

	int currentElementPosition = 0;
	for ( int iSection = 0; iSection < nTotalSections; ++ iSection )
	{
		CgnsSection & cgnsSection = cgnsZsection.GetCgnsSection( iSection );
		SectionMarker & section = iSection < nVolSections
			? *volSec.data[ iSection ]
			: *bcSec.data[ iSection - nVolSections ];

		int nElements = section.nElements;
		cgnsSection.SetSectionInfo( section.name, section.cgns_type, currentElementPosition + 1, currentElementPosition + nElements );
		cgnsSection.CreateConnList();
		currentElementPosition += nElements;

		int position = 0;
		for ( int iElement = 0; iElement < nElements; ++ iElement )
		{
			IntField & element = section.elements[ iElement ];
			int nElementNodes = element.size();
			for ( int iNode = 0; iNode < nElementNodes; ++ iNode )
			{
				cgnsSection.connList[ position ++ ]= element[ iNode ] + 1;
			}
		}
	}

	for ( int iSection = 0; iSection < nTotalSections; ++ iSection )
	{
		CgnsSection & cgnsSection = cgnsZsection.GetCgnsSection( iSection );
		cgnsSection.SetElemPosition();
	}

	CgnsZbc & cgnsZbc = cgnsZone.RequireCgnsZbc();
	CgnsZbcBoco & cgnsZbcBoco = cgnsZbc.RequireCgnsZbcBoco();
	cgnsZbcBoco.ReadZnboco( scalarBccos->bccos.size() );
	cgnsZbcBoco.CreateCgnsZbc();

	int currentBcElementPosition = nVolCell;
	for ( int iBcco = 0; iBcco < cgnsZbcBoco.GetNBoco(); ++ iBcco )
	{
		CgnsBcBoco & cgnsBcBoco = cgnsZbcBoco.GetCgnsBc( iBcco );
		ScalarBcco & scalarBcco = *scalarBccos->bccos[ iBcco ];
		int nElements = scalarBcco.eTypes.GetNElements();
		cgnsBcBoco.name = scalarBcco.bcName;
		cgnsBcBoco.nElements = scalarBcco.eTypes.GetNElements();
		cgnsBcBoco.bcType = static_cast< BCType_t >( scalarBcco.bcType );
		cgnsBcBoco.pointSetType = PointList;
		cgnsBcBoco.CreateCgnsBcBoco();

		for ( int iElement = 0; iElement < nElements; ++ iElement )
		{
			int bcElemId = scalarBcco.local_globalIds[ iElement ];
			cgnsBcBoco.SetConnListValue( iElement, bcElemId + 1 + nVolCell );
		}
	}

}

void ScalarGrid::DumpCgnsGrid()
{
	std::string prjFileName = Prj::GetMeshOutputFileName( "scalar.cgns" );
	// Use stack allocation instead of raw pointer
	CgnsZbase cgnsZbase;
	cgnsZbase.nBases = 1;
	cgnsZbase.InitCgnsBase();

	for ( int iBase = 0; iBase < cgnsZbase.nBases; ++ iBase )
	{
		CgnsBase * cgnsBase = cgnsZbase.GetCgnsBase( iBase );

		cgnsBase->celldim = ONEFLOW::ONE_D;
		cgnsBase->phydim  = ONEFLOW::ONE_D;
		cgnsBase->baseName = ONEFLOW::AddString( "Base", cgnsBase->baseId );
		cgnsBase->nZones = 1;
		cgnsBase->AllocateAllCgnsZones();

		int iZone = 0;
		CgnsZone * cgnsZone = cgnsBase->GetCgnsZone( iZone );
		this->SetCgnsZone( *cgnsZone );
	}

	cgnsZbase.cgnsFile->OpenCgnsFile( prjFileName, CG_MODE_WRITE );
	cgnsZbase.DumpCgnsMultiBase();
	cgnsZbase.cgnsFile->CloseCgnsFile();
}

void ScalarGrid::GenerateGridFromCgns( const std::string & prjFileName )
{
	// Use stack allocation instead of raw pointer
	CgnsZbase cgnsZbase;
	cgnsZbase.OpenCgnsFile( prjFileName, CG_MODE_READ );
	cgnsZbase.ReadCgnsMultiBase();
	cgnsZbase.CloseCgnsFile();
	this->ReadFromCgnsZbase( cgnsZbase );
}

void ScalarGrid::ReadFromCgnsZbase( CgnsZbase & cgnsZbase )
{
	int iBase = 0;
	int iZone = 0;
	CgnsBase * cgnsBase = cgnsZbase.GetCgnsBase( 0 );
	CgnsZone * cgnsZone = cgnsBase->GetCgnsZone( iZone );
	this->ReadFromCgnsZone( *cgnsZone );
}

void ScalarGrid::ReadFromCgnsZone( CgnsZone & cgnsZone )
{
	// Stage the complete import so malformed input cannot destroy the current mesh.
	ScalarGrid importedGrid;

	// Load coordinates first because PushElement validates one-based CGNS node IDs
	// against the imported mesh's current node count.
	CgnsCoor & cgnsCoor = cgnsZone.RequireCgnsCoor();
	NodeMesh & nodeMesh = cgnsCoor.RequireNodeMesh();
	if ( nodeMesh.xN.size() != nodeMesh.yN.size() ||
		 nodeMesh.xN.size() != nodeMesh.zN.size() )
	{
		throw std::invalid_argument( "ScalarGrid::ReadFromCgnsZone: coordinate arrays have inconsistent sizes" );
	}

	for ( size_t i = 0; i < nodeMesh.xN.size(); ++ i )
	{
		importedGrid.xn.AddData( nodeMesh.xN[ i ] );
		importedGrid.yn.AddData( nodeMesh.yN[ i ] );
		importedGrid.zn.AddData( nodeMesh.zN[ i ] );
	}


	std::cout << "   Convert Cgns Section Data to ScalarGrid......\n";
	std::cout << "\n";
	CgnsZsection & cgnsZsection = cgnsZone.RequireCgnsZsection();
	const int nSections = cgnsZsection.GetNSections();
	for ( int iSection = 0; iSection < nSections; ++ iSection )
	{
		std::cout << "-->iSection     = " << iSection << " numberOfCgnsSections = " << nSections << "\n";
		CgnsSection & cgnsSection = cgnsZsection.GetCgnsSection( iSection );

		const bool isHomogeneousVolumeSection =
			ONEFLOW::IsBasicVolumeElementType( cgnsSection.eType );
		const bool isMixedSection = cgnsSection.eType == MIXED;
		if ( ! isHomogeneousVolumeSection && ! isMixedSection ) continue;

		if ( cgnsSection.nElement < 0 ||
			 cgnsSection.eTypeList.size() != static_cast< size_t >( cgnsSection.nElement ) ||
			 cgnsSection.ePosList.size() != static_cast< size_t >( cgnsSection.nElement ) + 1 )
		{
			throw std::invalid_argument(
				"ScalarGrid::ReadFromCgnsZone: inconsistent element metadata in volume section" );
		}

		for ( int iElem = 0; iElem < cgnsSection.nElement; ++ iElem )
		{
			const int eType = cgnsSection.eTypeList[ iElem ];
			if ( ! ONEFLOW::IsBasicVolumeElementType( eType ) ) continue;

			const int nodeCount = ONEFLOW::GetElementNodeNumbers( eType );
			const long long connectionBegin = static_cast< long long >( cgnsSection.ePosList[ iElem ] ) +
				( isMixedSection ? 1LL : static_cast< long long >( cgnsSection.pos_shift ) );
			const long long nextConnectionBegin = static_cast< long long >( cgnsSection.ePosList[ iElem + 1 ] ) +
				( isMixedSection ? 0LL : static_cast< long long >( cgnsSection.pos_shift ) );
			if ( nodeCount <= 0 ||
				 connectionBegin < 0 ||
				 nextConnectionBegin < connectionBegin ||
				 nextConnectionBegin > static_cast< long long >( cgnsSection.connList.size() ) ||
				 connectionBegin + nodeCount > nextConnectionBegin )
			{
				throw std::invalid_argument(
					"ScalarGrid::ReadFromCgnsZone: invalid connectivity range in volume section" );
			}

			CgIntField eNodeId;
			cgnsSection.GetElementNodeId( iElem, eNodeId );
			importedGrid.PushElement( eNodeId, eType );
		}
	}


	const size_t nodeCount = importedGrid.xn.GetNElements();
	for ( const std::vector< int > & element : importedGrid.elements.data )
	{
		for ( const int nodeId : element )
		{
			if ( nodeId < 0 || static_cast< size_t >( nodeId ) >= nodeCount )
			{
				throw std::invalid_argument( "ScalarGrid::ReadFromCgnsZone: connectivity references a node outside the imported zone" );
			}
		}
	}

	// Commit only after sections, coordinates, and connectivity all validate.
	this->ResetMeshData();
	this->elements.data = std::move( importedGrid.elements.data );
	this->eTypes.data = std::move( importedGrid.eTypes.data );
	this->xn.data = std::move( importedGrid.xn.data );
	this->yn.data = std::move( importedGrid.yn.data );
	this->zn.data = std::move( importedGrid.zn.data );

	// Keep cached counts consistent with the newly committed storage.
	this->nNodes = this->GetNNodes();
	this->nCells = this->GetNCells();
}

void ScalarGrid::PushElement( CgIntField & eNodeId, int eType )
{
	if ( eType < 0 || eType >= ElementHome::GetUnitElementSlotCount() ||
		 ! ONEFLOW::IsBasicVolumeElementType( eType ) )
	{
		throw std::invalid_argument( "ScalarGrid::PushElement: invalid volume element type" );
	}

	const int expectedNodeCount = ONEFLOW::GetElementNodeNumbers( eType );
	if ( expectedNodeCount <= 0 || eNodeId.size() != static_cast< size_t >( expectedNodeCount ) )
	{
		throw std::invalid_argument( "ScalarGrid::PushElement: connectivity does not match the element type" );
	}

	// CGNS connectivity is one-based; validate it before converting to internal indices.
	const int nodeCount = this->GetNNodes();
	for ( int nodeId : eNodeId )
	{
		if ( nodeId <= 0 || nodeId > nodeCount )
		{
			throw std::invalid_argument( "ScalarGrid::PushElement: CGNS node index is outside the node range" );
		}
	}

	// Cell connectivity defines the topology's source data; discard derived state
	// before appending a cell to an existing mesh.
	this->ResetTopologyData();

	IntList elem;
	for ( int nodeId : eNodeId )
	{
		elem.AddData( nodeId - 1 );
	}

	this->elements.AddElem( elem );
	this->eTypes.AddData( eType );
}

void ScalarGrid::PushElement( int p1, int p2, int eType )
{
	if ( eType < 0 || eType >= ElementHome::GetUnitElementSlotCount() ||
		 ONEFLOW::GetElementNodeNumbers( eType ) != 2 )
	{
		throw std::invalid_argument( "ScalarGrid::PushElement: element type must describe a two-node element" );
	}

	const int nodeCount = this->GetNNodes();
	if ( p1 < 0 || p1 >= nodeCount || p2 < 0 || p2 >= nodeCount )
	{
		throw std::out_of_range( "ScalarGrid::PushElement: node index is outside the mesh node range" );
	}
	if ( p1 == p2 )
	{
		throw std::invalid_argument( "ScalarGrid::PushElement: cell connectivity must reference distinct nodes" );
	}

	// Validate the source connectivity before invalidating derived topology.
	this->ResetTopologyData();

	IntList elem;
	elem.AddData( p1 );
	elem.AddData( p2 );

	this->elements.AddElem( elem );
	this->eTypes.AddData( eType );
}

void ScalarGrid::PushBoundaryFace( int pt, int eType )
{
	this->ResetTopologyData();

	IntList elem;
	elem.AddData( pt );

	this->boundaryElements.AddElem( elem );
	this->bcETypes.AddData( eType );
}

void ScalarGrid::AllocGeom()
{
	this->nFaces = this->GetNFaces();
	this->nCells = this->GetNCells();
	this->nBFaces = this->GetNBFaces();

	this->nTCells = this->GetNTCells();

	this->xfc.Resize( this->nFaces );
	this->yfc.Resize( this->nFaces );
	this->zfc.Resize( this->nFaces );

	this->xfn.Resize( this->nFaces );
	this->yfn.Resize( this->nFaces );
	this->zfn.Resize( this->nFaces );
	this->area.Resize( this->nFaces );

	this->xcc.Resize( nTCells );
	this->ycc.Resize( nTCells );
	this->zcc.Resize( nTCells );
	this->vol.Resize( nTCells );
}

void ScalarGrid::CalcMetrics1D()
{
	this->ResetGeometryData();

	//must compute face center first for one dimensional case
	//then face normal
	this->AllocGeom();
	this->CalcFaceCenter1D();
	this->CalcCellCenterVol1D();
	this->CalcFaceNormal1D();
	this->CalcGhostCellCenterVol1D();
}

void ScalarGrid::CalcFaceCenter1D()
{
	this->nFaces = this->GetNFaces();
	for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
	{
		std::vector< int > & faceNodes = this->faces[ iFace ];
		int p1 = faceNodes[ 0 ];
		int p2 = faceNodes[ 0 ];
		this->xfc[ iFace ] = half * ( this->xn[ p1 ] + this->xn[ p2 ] );
		this->yfc[ iFace ] = half * ( this->yn[ p1 ] + this->yn[ p2 ] );
		this->zfc[ iFace ] = half * ( this->zn[ p1 ] + this->zn[ p2 ] );
	}
}

//void ScalarGrid::CalcCellCenterVol1D()
//{
//	this->nCells = this->GetNCells();
//	
//	for ( size_t iCell = 0; iCell < this->nCells; ++ iCell )
//	{
//		std::vector< int > & element = this->elements[ iCell ];
//		int p1 = element[ 0 ];
//		int p2 = element[ 1 ];
//		this->xcc[ iCell  ] = half * ( this->xn[ p1 ] + this->xn[ p2 ] );
//		this->ycc[ iCell  ] = half * ( this->yn[ p1 ] + this->yn[ p2 ] );
//		this->zcc[ iCell  ] = half * ( this->zn[ p1 ] + this->zn[ p2 ] );
//		Real dx = this->xn[ p2 ] - this->xn[ p1 ];
//		Real dy = this->yn[ p2 ] - this->yn[ p1 ];
//		Real dz = this->zn[ p2 ] - this->zn[ p1 ];
//		this->vol[ iCell  ] = ONEFLOW::DIST( dx, dy, dz );
//	}
//}


void ScalarGrid::CalcCellCenter1D()
{
	this->xcc = 0;
	this->ycc = 0;
	this->zcc = 0;

	this->nFaces = this->GetNFaces();
	this->nBFaces = this->GetNBFaces();
	for ( int iFace = 0; iFace < this->nBFaces; ++ iFace )
	{
		std::vector< int > & faceNodes = this->faces[ iFace ];
		int lc  = this->lc[ iFace ];

		int pt = faceNodes[ 0 ];

		this->xcc[ lc ] += this->xn[ pt ];
		this->ycc[ lc ] += this->yn[ pt ];
		this->zcc[ lc ] += this->zn[ pt ];
	}

	for ( int iFace = nBFaces; iFace < this->nFaces; ++ iFace )
	{
		std::vector< int > & faceNodes = this->faces[ iFace ];
		int lc  = this->lc[ iFace ];
		int rc  = this->rc[ iFace ];

		int pt = faceNodes[ 0 ];

		this->xcc[ lc ] += this->xn[ pt ];
		this->ycc[ lc ] += this->yn[ pt ];
		this->zcc[ lc ] += this->zn[ pt ];

		this->xcc[ rc ] += this->xn[ pt ];
		this->ycc[ rc ] += this->yn[ pt ];
		this->zcc[ rc ] += this->zn[ pt ];
	}

	this->nCells = this->GetNCells();

	for ( size_t iCell = 0; iCell < this->nCells; ++ iCell )
	{
		this->xcc[ iCell  ] *= half;
		this->ycc[ iCell  ] *= half;
		this->zcc[ iCell  ] *= half;
	}
}

void ScalarGrid::CalcCellVolume1D()
{
	this->vol = 0;

	this->nFaces = this->GetNFaces();
	for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
	{
		std::vector< int > & faceNodes = this->faces[ iFace ];
		int lc  = this->lc[ iFace ];
		int rc  = this->rc[ iFace ];

		int pt = faceNodes[ 0 ];

		Real dxl = this->xfc[ iFace ] - this->xcc[ lc ];
		Real dyl = this->yfc[ iFace ] - this->ycc[ lc ];
		Real dzl = this->zfc[ iFace ] - this->zcc[ lc ];

		Real dxr = this->xfc[ iFace ] - this->xcc[ rc ];
		Real dyr = this->yfc[ iFace ] - this->ycc[ rc ];
		Real dzr = this->zfc[ iFace ] - this->zcc[ rc ];

		Real dsl = ONEFLOW::DIST( dxl, dyl, dzl );
		Real dsr = ONEFLOW::DIST( dxr, dyr, dzr );

		this->vol[ lc ] += dsl;
		this->vol[ rc ] += dsr;
	}
}

void ScalarGrid::CalcCellCenterVol1D()
{
	this->CalcCellCenter1D();
	this->CalcCellVolume1D();
}

void ScalarGrid::CalcFaceNormal1D()
{
	this->nFaces = this->GetNFaces();
	for ( size_t iFace = 0; iFace < this->nFaces; ++ iFace )
	{
		int lc  = this->lc[ iFace ];

		Real dx = this->xfc[ iFace ] - this->xcc[ lc ];
		Real dy = this->yfc[ iFace ] - this->ycc[ lc ];
		Real dz = this->zfc[ iFace ] - this->zcc[ lc ];
		Real ds = ONEFLOW::DIST( dx, dy, dz );

		Real factor   = 1.0 / ( ds + SMALL );
		this->xfn[ iFace ] = factor * dx;
		this->yfn[ iFace ] = factor * dy;
		this->zfn[ iFace ] = factor * dz;

		this->area[ iFace ] = 1.0;
	}
}

void ScalarGrid::CalcGhostCellCenterVol1D()
{
	this->nBFaces = this->GetNBFaces();
	this->nCells = this->GetNCells();

	// For ghost cells
	for ( size_t iFace = 0; iFace < nBFaces; ++ iFace )
	{
		int lc = this->lc[ iFace ];
		int rc = this->rc[ iFace ];
		if ( this->area[ iFace ] > SMALL )
		{
			Real tmp = 2.0 * ( ( this->xcc[ lc ] - this->xfc[ iFace ] ) * this->xfn[ iFace ]
				             + ( this->ycc[ lc ] - this->yfc[ iFace ] ) * this->yfn[ iFace ]
				             + ( this->zcc[ lc ] - this->zfc[ iFace ] ) * this->zfn[ iFace ] );
			this->xcc[ rc ] = this->xcc[ lc ] - this->xfn[ iFace ] * tmp;
			this->ycc[ rc ] = this->ycc[ lc ] - this->yfn[ iFace ] * tmp;
			this->zcc[ rc ] = this->zcc[ lc ] - this->zfn[ iFace ] * tmp;
		}
		else
		{
			// Degenerated faces
			this->xcc[ rc ] = - this->xcc[ lc ] + 2.0 * this->xfc[ iFace ];
			this->ycc[ rc ] = - this->ycc[ lc ] + 2.0 * this->yfc[ iFace ];
			this->zcc[ rc ] = - this->zcc[ lc ] + 2.0 * this->zfc[ iFace ];
		}
		this->vol[ rc ] = this->vol[ lc ];
	}
}

std::vector< IntSet > ScalarGrid::CollectBoundaryVertexSets( int nodeCount ) const
{
	if ( this->scalarBccos == nullptr )
	{
		throw std::logic_error( "ScalarGrid::CalcTopology: boundary condition collection is not initialized" );
	}

	std::vector< IntSet > boundaryVertexSets;
	boundaryVertexSets.reserve( this->scalarBccos->bccos.size() );
	for ( const std::unique_ptr< ScalarBcco > & boundaryCondition : this->scalarBccos->bccos )
	{
		if ( boundaryCondition == nullptr )
		{
			throw std::runtime_error( "ScalarGrid::CalcTopology: boundary condition collection contains a null entry" );
		}

		IntSet boundaryVertices;
		for ( int nodeId : boundaryCondition->vertexList )
		{
			if ( nodeId < 0 || nodeId >= nodeCount )
			{
				throw std::runtime_error( "ScalarGrid::CalcTopology: boundary condition references an invalid node index" );
			}
			boundaryVertices.insert( nodeId );
		}
		boundaryVertexSets.push_back( std::move( boundaryVertices ) );
	}
	return boundaryVertexSets;
}

void ScalarGrid::CalcTopology()
{
	const int nodeCount = this->GetNNodes();
	const int cellCount = this->GetNCells();

	if ( this->yn.GetNElements() != static_cast< size_t >( nodeCount ) ||
		 this->zn.GetNElements() != static_cast< size_t >( nodeCount ) )
	{
		throw std::runtime_error( "ScalarGrid::CalcTopology: coordinate arrays have inconsistent sizes" );
	}

	if ( this->elements.GetNElements() != static_cast< size_t >( cellCount ) )
	{
		throw std::runtime_error( "ScalarGrid::CalcTopology: cell connectivity and element type arrays have inconsistent sizes" );
	}

	for ( int iCell = 0; iCell < cellCount; ++ iCell )
	{
		const int elementType = this->eTypes[ iCell ];
		if ( elementType < 0 || elementType >= ElementHome::GetUnitElementSlotCount() )
		{
			throw std::runtime_error( "ScalarGrid::CalcTopology: invalid element type" );
		}

		const int expectedNodeCount = ONEFLOW::GetElementNodeNumbers( elementType );
		const std::vector< int > & element = this->elements[ iCell ];
		if ( expectedNodeCount <= 0 || element.size() != static_cast< size_t >( expectedNodeCount ) )
		{
			throw std::runtime_error( "ScalarGrid::CalcTopology: cell connectivity does not match its element type" );
		}

		for ( int nodeId : element )
		{
			if ( nodeId < 0 || nodeId >= nodeCount )
			{
				throw std::runtime_error( "ScalarGrid::CalcTopology: cell references an invalid node index" );
			}
		}

		std::vector< int > sortedNodeIds( element );
		std::sort( sortedNodeIds.begin(), sortedNodeIds.end() );
		if ( std::adjacent_find( sortedNodeIds.begin(), sortedNodeIds.end() ) != sortedNodeIds.end() )
		{
			throw std::runtime_error( "ScalarGrid::CalcTopology: cell connectivity contains duplicate node indices" );
		}
	}

	const std::vector< IntSet > boundaryVertexSets = this->CollectBoundaryVertexSets( nodeCount );

	// A face in a conforming volume/line mesh may belong to at most two cells.
	// Check incidence before resetting the current topology, so malformed meshes
	// cannot silently overwrite a third cell or destroy the previous topology.
	HXLookup< int > incidenceLookup;
	std::vector< int > faceOwnerCell;
	std::vector< int > faceIncidenceCount;
	std::vector< std::vector< int > > faceNodeLists;
	for ( int iCell = 0; iCell < cellCount; ++ iCell )
	{
		const std::vector< int > & element = this->elements[ iCell ];
		const UnitElement & elementDefinition = ElementHome::GetUnitElementDefinition( this->eTypes[ iCell ] );
		const int localFaceCount = elementDefinition.GetElementFaceCount();
		for ( int iLocalFace = 0; iLocalFace < localFaceCount; ++ iLocalFace )
		{
			const IntField & localFaceNodes = elementDefinition.GetElementFace( iLocalFace );
			IntField faceNodes;
			faceNodes.reserve( localFaceNodes.size() );
			for ( int localNodeId : localFaceNodes )
			{
				faceNodes.push_back( element[ localNodeId ] );
			}

			auto [ faceIndex, isNew ] = incidenceLookup.FindOrAdd( faceNodes );
			if ( isNew )
			{
				faceOwnerCell.push_back( iCell );
				faceIncidenceCount.push_back( 1 );
				faceNodeLists.push_back( std::move( faceNodes ) );
				continue;
			}

			if ( faceOwnerCell[ faceIndex ] == iCell )
			{
				throw std::runtime_error( "ScalarGrid::CalcTopology: a cell contains a duplicate face" );
			}
			if ( faceIncidenceCount[ faceIndex ] >= 2 )
			{
				throw std::runtime_error( "ScalarGrid::CalcTopology: a face is shared by more than two cells" );
			}
			++ faceIncidenceCount[ faceIndex ];
		}
	}

	// When boundary-condition metadata is supplied, require each exterior
	// face to match exactly one condition. Some import paths construct scalar
	// topology before importing BC metadata, so an empty collection is valid
	// at this stage and must not prevent topology construction.
	if ( ! boundaryVertexSets.empty() )
	{
		for ( size_t iFace = 0; iFace < faceIncidenceCount.size(); ++ iFace )
		{
			if ( faceIncidenceCount[ iFace ] != 1 )
			{
				continue;
			}

			int matchingBoundaryConditions = 0;
			for ( const IntSet & boundaryVertices : boundaryVertexSets )
			{
				if ( this->CheckBcFace( boundaryVertices, faceNodeLists[ iFace ] ) )
				{
					++ matchingBoundaryConditions;
					if ( matchingBoundaryConditions > 1 )
					{
						throw std::runtime_error( "ScalarGrid::CalcTopology: boundary face matches multiple boundary conditions" );
					}
				}
			}
			if ( matchingBoundaryConditions == 0 )
			{
				throw std::runtime_error( "ScalarGrid::CalcTopology: boundary face does not match any boundary condition" );
			}
		}
	}

	// Validate the input before clearing existing topology so failed rebuilds
	// do not destroy a previously available topology.
	this->ResetTopologyData();

	this->nNodes = nodeCount;
	this->nCells = cellCount;

	// Use HXLookup to manage unique faces (key is automatically sorted)
	HXLookup<int> faceLookup;

	// Reserve for all element-local faces to avoid repeated growth.
	int estimatedFaces = 0;
	for ( int iCell = 0; iCell < nCells; ++ iCell )
	{
		int eType = this->eTypes[ iCell ];
		const UnitElement & elementDefinition = ElementHome::GetUnitElementDefinition( eType );
		estimatedFaces += elementDefinition.GetElementFaceCount();
	}

	this->lc.Reserve(estimatedFaces);
	this->rc.Reserve(estimatedFaces);
	this->lpos.Reserve(estimatedFaces);
	this->rpos.Reserve(estimatedFaces);
	this->fTypes.Reserve(estimatedFaces);
	this->fBcTypes.Reserve(estimatedFaces);
	this->faces.Reserve(estimatedFaces);

	for ( int iCell = 0; iCell < nCells; ++ iCell )
	{
		const std::vector<int>& element = elements[iCell];
		int eType = eTypes[iCell];
		const UnitElement & elementDefinition = ElementHome::GetUnitElementDefinition(eType);
		const int numberOfFaceInElement = elementDefinition.GetElementFaceCount();

		for ( int iLocalFace = 0; iLocalFace < numberOfFaceInElement; ++ iLocalFace )
		{
			const IntField& localFaceNodeIndexArray = elementDefinition.GetElementFace(iLocalFace);
			const int faceType = elementDefinition.GetFaceType(iLocalFace);

			// Build the global node array for the current face
			IntField faceNodeIndexArray;
			faceNodeIndexArray.reserve(localFaceNodeIndexArray.size());
			for (int nodeIndex : localFaceNodeIndexArray)
			{
				faceNodeIndexArray.push_back(element[nodeIndex]);
			}

			// Find or add the face (HXLookup automatically sorts the nodes)
			auto [faceIndex, isNew] = faceLookup.FindOrAdd(faceNodeIndexArray);

			if ( isNew )
			{
				// New face: add face data
				this->lc.data.push_back(iCell);
				this->rc.data.push_back(ONEFLOW::INVALID_INDEX);
				this->lpos.data.push_back(iLocalFace);
				this->rpos.data.push_back(ONEFLOW::INVALID_INDEX);
				this->fTypes.data.push_back(faceType);
				this->fBcTypes.data.push_back(ONEFLOW::INVALID_INDEX);
				this->faces.data.push_back(std::move(faceNodeIndexArray));
			}
			else
			{
				// Existing face: update the right cell information
				this->fBcTypes[faceIndex] = ONEFLOW::BCTypeNull;  // Internal face
				this->rc[faceIndex] = iCell;
				this->rpos[faceIndex] = iLocalFace;
			}
		}
	}

	this->ReorderFaces();
	this->ScanBcFace();
	this->SetBcGhostCell();
}

void ScalarGrid::ScanBcFace()
{
	this->AllocateBc();
	scalarBccos->ScanBcFace( *this );
	this->SetBcTypes();
}

void ScalarGrid::CalcOrderMap( IntList &orderMap )
{
	this->nFaces = this->GetNFaces();
	orderMap.Resize( this->nFaces );

	int iBoundaryFaceCount = 0;
	int iCount = 0;
	for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
	{
		int rc = this->rc[ iFace ];
		if ( rc == ONEFLOW::INVALID_INDEX )
		{
			orderMap[ iCount ++ ] = iFace;
			++ iBoundaryFaceCount;
		}
	}

	this->nBFaces = iBoundaryFaceCount;

	for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
	{
		int rc = this->rc[ iFace ];
		if ( rc != ONEFLOW::INVALID_INDEX )
		{
			orderMap[ iCount ++ ] = iFace;
		}
	}

}

void ScalarGrid::ReorderFaces()
{
	IntList orderMap;
	this->CalcOrderMap( orderMap );

	const size_t nFaces = orderMap.GetNElements();
	if ( this->lc.GetNElements() != nFaces ||
		 this->rc.GetNElements() != nFaces ||
		 this->lpos.GetNElements() != nFaces ||
		 this->rpos.GetNElements() != nFaces ||
		 this->fTypes.GetNElements() != nFaces ||
		 this->fBcTypes.GetNElements() != nFaces ||
		 this->faces.GetNElements() != nFaces )
	{
		throw std::logic_error( "ScalarGrid::ReorderFaces: face topology arrays have inconsistent sizes" );
	}

	this->lc.ReOrder( orderMap );
	this->rc.ReOrder( orderMap );
	this->lpos.ReOrder( orderMap );
	this->rpos.ReOrder( orderMap );
	this->fTypes.ReOrder( orderMap );
	this->fBcTypes.ReOrder( orderMap );
	this->faces.ReOrder( orderMap );
}

void ScalarGrid::SetBcGhostCell()
{
	this->nCells = this->GetNCells();
	this->nBFaces = this->GetNBFaces();

	for ( int iFace = 0; iFace < this->nBFaces; ++ iFace )
	{
		this->rc[ iFace ] = iFace + nCells;
	}
}

bool ScalarGrid::CheckBcFace( const IntSet & bcVertex, const std::vector< int > & nodeId ) const
{
	int size = nodeId.size();
	for ( int iNode = 0; iNode < size; ++ iNode )
	{
		IntSet::iterator iter = bcVertex.find( nodeId[ iNode ] );
		if ( iter == bcVertex.end() )
		{
			return false;
		}
	}
	return true;
}

void ScalarGrid::AllocateBc()
{
	this->nFaces = this->GetNFaces();
	std::cout << " nFaces = " << nFaces << "\n";

	int nTraditionalBc = 0;
	for ( int iFace = 0; iFace < nFaces; ++ iFace )
	{
		int originalBcType = this->fBcTypes[ iFace ];
		if ( originalBcType == ONEFLOW::INVALID_INDEX )
		{
			++ nTraditionalBc;
		}
	}
	std::cout << " nTraditionalBc = " << nTraditionalBc << "\n";
	this->bcTypes.Resize( nTraditionalBc );
}

void ScalarGrid::ScanBcFace( IntSet& bcVertex, int bcType )
{
	this->nFaces = this->GetNFaces();
	int nBcFaces_local = 0;
	for ( int iFace = 0; iFace < nFaces; ++ iFace )
	{
		int originalBcType = this->fBcTypes[ iFace ];

		if ( originalBcType == ONEFLOW::INVALID_INDEX )
		{
			if ( this->CheckBcFace( bcVertex, this->faces[ iFace ] ) )
			{
				++ nBcFaces_local;

				this->fBcTypes[ iFace ] = bcType;
			}
		}
	}

	std::cout << " nFinalBcFace = " << nBcFaces_local << " bcType = " << bcType << std::endl;

}

void ScalarGrid::SetBcTypes()
{
	this->nBFaces = this->GetNBFaces();
	for ( int iFace = 0; iFace < nBFaces; ++ iFace )
	{
		int bcType = this->fBcTypes[ iFace ];
		this->bcTypes[ iFace ] = bcType;
	}
}

void ScalarGrid::CalcC2C( EList & c2c ) const
{
	// Rebuild adjacency from the current face connectivity on every call.
	const int nFaces = this->GetNFaces();
	const int nCells = this->GetNCells();
	const int nBFaces = this->GetNBFaces();
	if ( nBFaces < 0 || nBFaces > nFaces ||
		 this->lc.GetNElements() != static_cast< size_t >( nFaces ) ||
		 this->rc.GetNElements() != static_cast< size_t >( nFaces ) )
	{
		throw std::runtime_error( "ScalarGrid::CalcC2C: face topology arrays have inconsistent sizes" );
	}

	// Validate every cell index before using it to index an adjacency row.
	for ( int iFace = 0; iFace < nBFaces; ++ iFace )
	{
		if ( BC::IsInterfaceBc( this->bcTypes[ iFace ] ) )
		{
			const int leftCell = this->lc[ iFace ];
			if ( leftCell < 0 || leftCell >= nCells )
			{
				throw std::runtime_error( "ScalarGrid::CalcC2C: interface boundary face references an invalid physical cell" );
			}
		}
	}

	for ( int iFace = nBFaces; iFace < nFaces; ++ iFace )
	{
		const int leftCell = this->lc[ iFace ];
		const int rightCell = this->rc[ iFace ];
		if ( leftCell < 0 || leftCell >= nCells ||
			 rightCell < 0 || rightCell >= nCells || leftCell == rightCell )
		{
			throw std::runtime_error( "ScalarGrid::CalcC2C: internal face must connect two distinct valid physical cells" );
		}
	}

	// Build into a temporary so invalid topology never leaves partial adjacency.
	EList reconstructed;
	reconstructed.Resize( nCells );

	// If boundary is an INTERFACE, need to count ghost cell
	for ( int iFace = 0; iFace < nBFaces; ++ iFace )
	{
		if ( BC::IsInterfaceBc( this->bcTypes[ iFace ] ) )
		{
			reconstructed[ this->lc[ iFace ] ].push_back( this->rc[ iFace ] );
		}
	}

	for ( int iFace = nBFaces; iFace < nFaces; ++ iFace )
	{
		const int leftCell = this->lc[ iFace ];
		const int rightCell = this->rc[ iFace ];
		reconstructed[ leftCell ].push_back( rightCell );
		reconstructed[ rightCell ].push_back( leftCell );
	}

	c2c.data = std::move( reconstructed.data );
}

void ScalarGrid::CalcInterfaceToBcFace()
{
	const int nBFaces = this->GetNBFaces();
	const int nFaces = this->GetNFaces();
	if ( nBFaces > nFaces ||
		 this->lc.GetNElements() != static_cast< size_t >( nFaces ) ||
		 this->rc.GetNElements() != static_cast< size_t >( nFaces ) ||
		 this->fBcTypes.GetNElements() != static_cast< size_t >( nFaces ) )
	{
		throw std::runtime_error( "ScalarGrid::CalcInterfaceToBcFace: face topology arrays have inconsistent sizes" );
	}
	if ( this->scalarIFace == nullptr )
	{
		throw std::logic_error( "ScalarGrid::CalcInterfaceToBcFace: interface topology is not initialized" );
	}

	std::vector< int > interfaceToBcFace;
	if ( this->scalarIFace->GetNIFaces() == 0 )
	{
		this->scalarIFace->interface_to_bcface = std::move( interfaceToBcFace );
		return;
	}
	interfaceToBcFace.reserve( this->scalarIFace->GetNIFaces() );

	for ( int iBFace = 0; iBFace < nBFaces; ++ iBFace )
	{
		if ( ! BC::IsInterfaceBc( this->bcTypes[ iBFace ] ) )
		{
			continue;
		}

		interfaceToBcFace.push_back( iBFace );
	}

	if ( interfaceToBcFace.size() != static_cast< size_t >( this->scalarIFace->GetNIFaces() ) )
	{
		throw std::runtime_error( "ScalarGrid::CalcInterfaceToBcFace: interface face count does not match interface topology" );
	}

	// Commit the complete derived mapping only after validation succeeds.
	this->scalarIFace->interface_to_bcface = std::move( interfaceToBcFace );
}

void ScalarGrid::Normalize()
{
	const int nFaces = this->GetNFaces();
	const int nBFaces = this->GetNBFaces();
	const int nCells = this->GetNCells();
	if ( this->lc.GetNElements() != static_cast< size_t >( nFaces ) ||
		 this->rc.GetNElements() != static_cast< size_t >( nFaces ) ||
		 this->faces.GetNElements() != static_cast< size_t >( nFaces ) )
	{
		throw std::runtime_error( "ScalarGrid::Normalize: face topology arrays have inconsistent sizes" );
	}
	if ( nBFaces > nFaces )
	{
		throw std::runtime_error( "ScalarGrid::Normalize: boundary face count exceeds total face count" );
	}

	// Validate all physical adjacency before changing orientation or assigning ghost IDs.
	for ( int iFace = 0; iFace < nFaces; ++ iFace )
	{
		const int leftCell = this->lc[ iFace ];
		const int rightCell = this->rc[ iFace ];
		const bool leftIsPhysical = leftCell >= 0 && leftCell < nCells;
		const bool rightIsPhysical = rightCell >= 0 && rightCell < nCells;
		const bool leftIsMissing = leftCell == ONEFLOW::INVALID_INDEX;
		const bool rightIsMissing = rightCell == ONEFLOW::INVALID_INDEX;

		if ( iFace < nBFaces )
		{
			if ( leftIsPhysical == rightIsPhysical ||
				 ( ! leftIsPhysical && ! leftIsMissing ) ||
				 ( ! rightIsPhysical && ! rightIsMissing ) )
			{
				throw std::runtime_error(
					"ScalarGrid::Normalize: boundary face must have one physical cell and one missing neighbor" );
			}
		}
		else if ( ! leftIsPhysical || ! rightIsPhysical || leftCell == rightCell )
		{
			throw std::runtime_error(
				"ScalarGrid::Normalize: internal face must connect two distinct physical cells" );
		}
	}

	for ( int iFace = 0; iFace < nFaces; ++ iFace )
	{
		if ( this->lc[ iFace ] == ONEFLOW::INVALID_INDEX )
		{
			// Reverse face orientation so the physical cell is on the left.
			std::vector< int > & face = this->faces[ iFace ];
			std::reverse( face.begin(), face.end() );
			ONEFLOW::SWAP( this->lc[ iFace ], this->rc[ iFace ] );
		}
	}

	this->SetBcGhostCell();
}

void ScalarGrid::GetSId( int i_interface, int & sId ) const
{
	const auto & interfaceToBcFace = this->scalarIFace->interface_to_bcface;
	if ( i_interface < 0 || static_cast< size_t >( i_interface ) >= interfaceToBcFace.size() )
	{
		throw std::out_of_range( "ScalarGrid::GetSId: interface index is out of range" );
	}

	const int iBFace = interfaceToBcFace[ i_interface ];
	if ( iBFace < 0 || static_cast< size_t >( iBFace ) >= this->lc.GetNElements() )
	{
		throw std::out_of_range( "ScalarGrid::GetSId: boundary face index is out of range" );
	}
	sId = this->lc[ iBFace ];
}

void ScalarGrid::GetTId( int i_interface, int & tId )
{
	const auto & interfaceToBcFace = this->scalarIFace->interface_to_bcface;
	if ( i_interface < 0 || static_cast< size_t >( i_interface ) >= interfaceToBcFace.size() )
	{
		throw std::out_of_range( "ScalarGrid::GetTId: interface index is out of range" );
	}

	const int iBFace = interfaceToBcFace[ i_interface ];
	if ( iBFace < 0 || static_cast< size_t >( iBFace ) >= this->rc.GetNElements() )
	{
		throw std::out_of_range( "ScalarGrid::GetTId: boundary face index is out of range" );
	}
	tId = this->rc[ iBFace ];
}

void ScalarGrid::DumpCalcGrid()
{
	std::cout << "Dumping unstructured grid data files......\n";
	std::fstream file;
	const std::string fileName = Prj::GetMeshOutputFileName( "scalar.ofl" );
	Prj::CreateDirIfNeeded( fileName );
	Prj::OpenFile( file, fileName, std::ios_base::out | std::ios_base::binary );
	DataBook databook;
	this->WriteGrid( &databook );
	databook.WriteFile( file );
	Prj::CloseFile( file );
}

void ScalarGrid::WriteGrid( std::fstream & file )
{
	DataBook databook;
	this->WriteGrid( &databook );
	databook.WriteFile( file );
}

void ScalarGrid::WriteGrid( DataBook * databook )
{
	this->nNodes = this->GetNNodes();
	this->nCells = this->GetNCells();
	this->nFaces = this->GetNFaces();

	ONEFLOW::HXWrite( databook, this->nNodes );
	ONEFLOW::HXWrite( databook, this->nFaces );
	ONEFLOW::HXWrite( databook, this->nCells );

	std::cout << " number of nodes    : " << this->nNodes << std::endl;
	std::cout << " number of surfaces : " << this->nFaces << std::endl;
	std::cout << " number of elements : " << this->nCells << std::endl;

	//node
	ONEFLOW::HXWrite( databook, this->xn.data );
	ONEFLOW::HXWrite( databook, this->yn.data );
	ONEFLOW::HXWrite( databook, this->zn.data );

	std::cout << " dumping xn,yn,zn \n";

	ONEFLOW::HXWrite( databook, this->volBcType  );
	std::cout << " this->volBcType = " << this->volBcType << "\n";

	std::cout << " dumping eTypes \n";

	//element
	ONEFLOW::HXWrite( databook, this->eTypes.data );

	this->WriteGridFaceTopology( databook );
	this->WriteBoundaryTopology( databook );
}

void ScalarGrid::ReadCalcGrid()
{
	std::fstream file;
	const std::string fileName = Prj::GetMeshInputFileName( "scalar.ofl" );
	Prj::OpenFile( file, fileName, std::ios_base::in | std::ios_base::binary );
	DataBook databook;
	databook.ReadFile( file );
	this->ReadGrid( &databook );
	Prj::CloseFile( file );
}

void ScalarGrid::ReadGrid( std::fstream & file )
{
	DataBook databook;
	databook.ReadFile( file );
	this->ReadGrid( &databook );
}

void ScalarGrid::ReadGrid( DataBook * databook )
{
	if ( databook == nullptr )
	{
		throw std::invalid_argument( "ScalarGrid::ReadGrid: databook must not be null" );
	}

	std::cout << "Reading unstructured grid data files......\n";
	std::cout << "Grid dimension = " << Dim::dimension << std::endl;

	int nodeCount = 0;
	int faceCount = 0;
	int cellCount = 0;
	ONEFLOW::HXRead( databook, nodeCount );
	ONEFLOW::HXRead( databook, faceCount );
	ONEFLOW::HXRead( databook, cellCount );

	if ( nodeCount < 0 || faceCount < 0 || cellCount < 0 )
	{
		throw std::runtime_error( "ScalarGrid::ReadGrid: grid header contains a negative count" );
	}

	// Parse and validate the replacement independently; a failed read must not
	// leave the current grid with a mixture of old and partially loaded data.
	ScalarGrid stagedGrid;
	stagedGrid.CreateNodes( nodeCount );
	// CreateNodes invalidates topology-derived counts, so publish the file header
	// counts only after establishing the replacement node storage.
	stagedGrid.nNodes = nodeCount;
	stagedGrid.nFaces = faceCount;
	stagedGrid.nCells = cellCount;

	std::cout << " number of nodes    : " << stagedGrid.nNodes << std::endl;
	std::cout << " number of surfaces : " << stagedGrid.nFaces << std::endl;
	std::cout << " number of elements : " << stagedGrid.nCells << std::endl;
	std::cout << " Reading xn,yn,zn\n";
	ONEFLOW::HXRead( databook, stagedGrid.xn.data );
	ONEFLOW::HXRead( databook, stagedGrid.yn.data );
	ONEFLOW::HXRead( databook, stagedGrid.zn.data );

	std::cout << " Reading volBcType\n";
	stagedGrid.volBcType = -1000;
	ONEFLOW::HXRead( databook, stagedGrid.volBcType );
	std::cout << " this->volBcType = " << stagedGrid.volBcType << "\n";

	std::cout << " Reading eTypes\n";
	stagedGrid.eTypes.Resize( stagedGrid.nCells );
	ONEFLOW::HXRead( databook, stagedGrid.eTypes.data );

	std::cout << "The grid nodes have been read\n";
	stagedGrid.ReadGridFaceTopology( databook );
	stagedGrid.ReadBoundaryTopology( databook );
	stagedGrid.NormalizeBc();

	// Publish a complete, normalized mesh in one commit. The DataBase and grid
	// identity belong to this ScalarGrid instance and are intentionally retained.
	this->xn.data = std::move( stagedGrid.xn.data );
	this->yn.data = std::move( stagedGrid.yn.data );
	this->zn.data = std::move( stagedGrid.zn.data );
	this->xfc.data = std::move( stagedGrid.xfc.data );
	this->yfc.data = std::move( stagedGrid.yfc.data );
	this->zfc.data = std::move( stagedGrid.zfc.data );
	this->xfn.data = std::move( stagedGrid.xfn.data );
	this->yfn.data = std::move( stagedGrid.yfn.data );
	this->zfn.data = std::move( stagedGrid.zfn.data );
	this->area.data = std::move( stagedGrid.area.data );
	this->xcc.data = std::move( stagedGrid.xcc.data );
	this->ycc.data = std::move( stagedGrid.ycc.data );
	this->zcc.data = std::move( stagedGrid.zcc.data );
	this->vol.data = std::move( stagedGrid.vol.data );

	this->lc.data = std::move( stagedGrid.lc.data );
	this->rc.data = std::move( stagedGrid.rc.data );
	this->lpos.data = std::move( stagedGrid.lpos.data );
	this->rpos.data = std::move( stagedGrid.rpos.data );
	this->cell2faces = std::move( stagedGrid.cell2faces );
	this->c2fpos = std::move( stagedGrid.c2fpos );
	this->global_faceid = std::move( stagedGrid.global_faceid );
	this->faces.data = std::move( stagedGrid.faces.data );
	this->elements.data = std::move( stagedGrid.elements.data );
	this->boundaryElements.data = std::move( stagedGrid.boundaryElements.data );
	this->bcETypes.data = std::move( stagedGrid.bcETypes.data );
	this->fTypes.data = std::move( stagedGrid.fTypes.data );
	this->eTypes.data = std::move( stagedGrid.eTypes.data );
	this->fBcTypes.data = std::move( stagedGrid.fBcTypes.data );
	this->bcTypes.data = std::move( stagedGrid.bcTypes.data );
	this->bcNameIds.data = std::move( stagedGrid.bcNameIds.data );
	this->scalarBccos = std::move( stagedGrid.scalarBccos );
	this->scalarIFace = std::move( stagedGrid.scalarIFace );

	this->nNodes = stagedGrid.nNodes;
	this->nCells = stagedGrid.nCells;
	this->nBFaces = stagedGrid.nBFaces;
	this->nFaces = stagedGrid.nFaces;
	this->nTCells = stagedGrid.nTCells;
	this->volBcType = stagedGrid.volBcType;

	std::cout << "All the computing information is ready!\n";
}

void ScalarGrid::NormalizeBc()
{
	if ( this->nBFaces < 0 || this->nBFaces > this->nFaces ||
		 this->lc.GetNElements() != static_cast< size_t >( this->nFaces ) ||
		 this->rc.GetNElements() != static_cast< size_t >( this->nFaces ) )
	{
		throw std::logic_error( "ScalarGrid::NormalizeBc: inconsistent face topology sizes" );
	}

	// Validate the complete face adjacency before assigning any ghost-cell IDs.
	for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
	{
		const int leftCell = this->lc[ iFace ];
		const int rightCell = this->rc[ iFace ];

		if ( leftCell < 0 || leftCell >= this->nCells )
		{
			throw std::out_of_range( "ScalarGrid::NormalizeBc: invalid left cell index" );
		}

		if ( iFace < this->nBFaces )
		{
			const int expectedGhostCell = iFace + this->nCells;
			if ( rightCell != -1 && rightCell != expectedGhostCell )
			{
				throw std::logic_error(
					"ScalarGrid::NormalizeBc: boundary face must have one real adjacent cell" );
			}
		}
		else
		{
			if ( rightCell < 0 || rightCell >= this->nCells )
			{
				throw std::out_of_range( "ScalarGrid::NormalizeBc: invalid internal right cell index" );
			}
			if ( rightCell == leftCell )
			{
				throw std::logic_error(
					"ScalarGrid::NormalizeBc: internal face must connect distinct cells" );
			}
		}
	}

	for ( int iFace = 0; iFace < this->nBFaces; ++ iFace )
	{
		this->rc[ iFace ] = iFace + this->nCells;
	}
}

void ScalarGrid::CreateNodes( int numberOfNodes )
{
	if ( numberOfNodes < 0 )
	{
		throw std::invalid_argument( "ScalarGrid::CreateNodes: number of nodes must be non-negative" );
	}

	// Node storage defines the index space used by both topology and geometry.
	// Invalidate derived data before changing that index space.
	this->ResetTopologyData();

	this->xn.Resize( numberOfNodes );
	this->yn.Resize( numberOfNodes );
	this->zn.Resize( numberOfNodes );
}

void ScalarGrid::WriteGridFaceTopology( DataBook * databook )
{
	// Validate the complete topology before writing any part of the serialized record.
	this->ValidateFaceTopology( "ScalarGrid::WriteGridFaceTopology" );

	std::cout << " Dumping this->fTypes \n";
	ONEFLOW::HXWrite( databook, this->fTypes.data );

	//std::cout << "fTypes = \n";
	//for ( int iFace = 0; iFace < this->fTypes.data.size(); ++ iFace )
	//{
	//	std::cout << this->fTypes.data[ iFace ] << " ";
	//}
	//std::cout << "\n";

	IntField numFaceNode( this->nFaces );

	for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
	{
		numFaceNode[ iFace ] = this->faces[ iFace ].size();
	}

	std::cout << " Dumping numFaceNode \n";

	ONEFLOW::HXWrite( databook, numFaceNode );

	int nsum = ONEFLOW::SUM( numFaceNode );
	std::cout << " nsum = " << nsum << "\n";
	std::cout << "numFaceNode = \n";
	//for ( int iFace = 0; iFace < numFaceNode.size(); ++ iFace )
	//{
	//	std::cout << numFaceNode[ iFace ] << " ";
	//}
	//std::cout << "\n";

	IntField faceNodeMem;
	faceNodeMem.reserve( nsum );

	for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
	{
		int nNodes = numFaceNode[ iFace ];
		for ( int iNode = 0; iNode < nNodes; ++ iNode )
		{
			faceNodeMem.push_back( this->faces[ iFace ][ iNode ] );
		}
	}
	std::cout << " Dumping faceNodeMem \n";
	ONEFLOW::HXWrite( databook, faceNodeMem );

	ONEFLOW::HXWrite( databook, this->lc.data );

	// Ghost-cell IDs are reconstructed by NormalizeBc; persist boundary sentinels instead.
	IntField serializedRightCells = this->rc.data;
	for ( int iFace = 0; iFace < this->nBFaces; ++ iFace )
	{
		serializedRightCells[ iFace ] = -1;
	}
	ONEFLOW::HXWrite( databook, serializedRightCells );
}

void ScalarGrid::ReadGridFaceTopology( DataBook * databook )
{
	// Parse and validate into local storage so malformed input cannot leave partial topology.
	IntList faceTypes;
	faceTypes.Resize( this->nFaces );

	std::cout << " Reading this->fTypes\n";
	ONEFLOW::HXRead( databook, faceTypes.data );

	IntField numFaceNode( this->nFaces );
	std::cout << " Reading numFaceNode\n";
	ONEFLOW::HXRead( databook, numFaceNode );

	size_t nsum = 0;
	for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
	{
		const int nFaceNodes = numFaceNode[ iFace ];
		// Every face has at least one node, including point-faces in 1D.
		if ( nFaceNodes <= 0 ||
			 static_cast< size_t >( nFaceNodes ) > std::numeric_limits< int >::max() - nsum )
		{
			throw std::runtime_error( "ScalarGrid::ReadGridFaceTopology: invalid or overflowing face-node count" );
		}
		nsum += static_cast< size_t >( nFaceNodes );
	}

	std::cout << " nsum = " << nsum << "\n";
	std::cout << " this->nFaces = " << this->nFaces << "\n";
	IntField faceNodeMem( nsum );
	std::cout << " Reading faceNodeMem\n";
	ONEFLOW::HXRead( databook, faceNodeMem );

	for ( const int nodeId : faceNodeMem )
	{
		if ( nodeId < 0 || nodeId >= this->nNodes )
		{
			throw std::runtime_error( "ScalarGrid::ReadGridFaceTopology: face references an invalid node index" );
		}
	}

	int faceNodeOffset = 0;
	for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
	{
		std::unordered_set< int > uniqueNodeIds;
		const int nFaceNodes = numFaceNode[ iFace ];
		for ( int iNode = 0; iNode < nFaceNodes; ++ iNode )
		{
			const int nodeId = faceNodeMem[ faceNodeOffset + iNode ];
			if ( ! uniqueNodeIds.insert( nodeId ).second )
			{
				throw std::runtime_error( "ScalarGrid::ReadGridFaceTopology: face contains duplicate node indices" );
			}
		}
		faceNodeOffset += nFaceNodes;
	}

	EList stagedFaces;
	stagedFaces.Resize( this->nFaces );
	int ipos = 0;
	for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
	{
		std::vector< int > & face = stagedFaces[ iFace ];
		const int nNodes = numFaceNode[ iFace ];
		face.reserve( nNodes );
		for ( int iNode = 0; iNode < nNodes; ++ iNode )
		{
			face.push_back( faceNodeMem[ ipos ++ ] );
		}
	}

	IntList leftCells;
	IntList rightCells;
	leftCells.Resize( this->nFaces );
	rightCells.Resize( this->nFaces );

	std::cout << "Setting the connection mode of face to cell......\n";
	ONEFLOW::HXRead( databook, leftCells.data );
	ONEFLOW::HXRead( databook, rightCells.data );

	// Validate all cell references before changing face orientation.
	for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
	{
		const int leftCell = leftCells[ iFace ];
		const int rightCell = rightCells[ iFace ];
		if ( leftCell < -1 || rightCell < -1 ||
			 leftCell >= this->nCells || rightCell >= this->nCells ||
			 ( leftCell < 0 && rightCell < 0 ) )
		{
			throw std::runtime_error( "ScalarGrid::ReadGridFaceTopology: face references invalid adjacent cells" );
		}
		if ( leftCell >= 0 && leftCell == rightCell )
		{
			throw std::runtime_error( "ScalarGrid::ReadGridFaceTopology: face must not connect a cell to itself" );
		}
	}

	for ( int iFace = 0; iFace < this->nFaces; ++ iFace )
	{
		if ( leftCells[ iFace ] < 0 )
		{
			std::reverse( stagedFaces[ iFace ].begin(), stagedFaces[ iFace ].end() );
			ONEFLOW::SWAP( leftCells[ iFace ], rightCells[ iFace ] );
		}
	}

	// Publish the complete topology only after all input and adjacency checks pass.
	this->faces.data = std::move( stagedFaces.data );
	this->fTypes.data = std::move( faceTypes.data );
	this->lc.data = std::move( leftCells.data );
	this->rc.data = std::move( rightCells.data );
}

void ScalarGrid::WriteBoundaryTopology( DataBook * databook )
{
	int nBFaces = this->GetNBFaces();
	ONEFLOW::HXWrite( databook, nBFaces );

	ONEFLOW::HXWrite( databook, this->bcTypes.data );
	this->bcNameIds = this->bcTypes;
	ONEFLOW::HXWrite( databook, this->bcNameIds.data );

	this->scalarIFace->WriteInterfaceTopology( databook );
}

void ScalarGrid::ReadBoundaryTopology( DataBook * databook )
{
	std::cout << "Setting the boundary condition......\n";
	int boundaryFaceCount = 0;
	ONEFLOW::HXRead( databook, boundaryFaceCount );

	if ( boundaryFaceCount < 0 || boundaryFaceCount > this->nFaces )
	{
		throw std::runtime_error( "ScalarGrid::ReadBoundaryTopology: boundary face count is outside the total face range" );
	}

	// Keep serialized boundary state local until the interface topology is consistent.
	IntList boundaryTypes;
	IntList boundaryNameIds;
	boundaryTypes.Resize( boundaryFaceCount );
	boundaryNameIds.Resize( boundaryFaceCount );
	ONEFLOW::HXRead( databook, boundaryTypes.data );
	ONEFLOW::HXRead( databook, boundaryNameIds.data );

	ScalarIFace stagedInterface;
	stagedInterface.ReadInterfaceTopology( databook );

	int nInterfaceBcFaces = 0;
	for ( int iBFace = 0; iBFace < boundaryFaceCount; ++ iBFace )
	{
		if ( BC::IsInterfaceBc( boundaryTypes[ iBFace ] ) )
		{
			++ nInterfaceBcFaces;
		}
	}
	if ( nInterfaceBcFaces != stagedInterface.GetNIFaces() )
	{
		throw std::runtime_error(
			"ScalarGrid::ReadBoundaryTopology: interface boundary-face count does not match interface topology" );
	}

	// Preserve runtime send/receive storage while replacing only serialized topology.
	if ( this->scalarIFace == nullptr )
		this->scalarIFace = std::make_unique< ScalarIFace >();
	this->scalarIFace->zones = std::move( stagedInterface.zones );
	this->scalarIFace->target_interfaces = std::move( stagedInterface.target_interfaces );
	this->scalarIFace->interface_to_bcface = std::move( stagedInterface.interface_to_bcface );
	this->scalarIFace->data = std::move( stagedInterface.data );
	this->scalarIFace->iglobalfaces.clear();
	this->scalarIFace->cells.clear();
	this->scalarIFace->global_to_local_interfaces.clear();
	this->scalarIFace->local_to_global_interfaces.clear();

	this->bcTypes.data = std::move( boundaryTypes.data );
	this->bcNameIds.data = std::move( boundaryNameIds.data );
	this->nBFaces = boundaryFaceCount;
}

//for partition
void ScalarGrid::AddPhysicalBcFace( int global_face_id, int bctype, int lcell, int rcell, int ftype )
{
	// Reserve every destination before changing the logical face record.
	ReserveForAppend( this->global_faceid );
	ReserveForAppend( this->bcTypes.data );
	ReserveForAppend( this->fBcTypes.data );
	ReserveForAppend( this->lc.data );
	ReserveForAppend( this->rc.data );
	ReserveForAppend( this->fTypes.data );

	this->global_faceid.push_back( global_face_id );
	this->bcTypes.AddData( bctype );
	this->fBcTypes.AddData( bctype );
	this->lc.AddData( lcell );
	this->rc.AddData( rcell );
	this->fTypes.AddData( ftype );
}

void ScalarGrid::AddInterfaceBcFace( int global_face_id, int bctype, int lcell, int rcell, int nei_zoneid, int nei_cellid, int ftype )
{
	// Reserve every destination before mutating interface or face topology.
	ReserveForAppend( this->global_faceid );
	ReserveForAppend( this->bcTypes.data );
	ReserveForAppend( this->fBcTypes.data );
	ReserveForAppend( this->lc.data );
	ReserveForAppend( this->rc.data );
	ReserveForAppend( this->fTypes.data );

	this->AddInterface( global_face_id, nei_zoneid, nei_cellid );

	this->global_faceid.push_back( global_face_id );
	this->bcTypes.AddData( bctype );
	this->fBcTypes.AddData( bctype );
	this->lc.AddData( lcell );
	this->rc.AddData( rcell );
	this->fTypes.AddData( ftype );
}

void ScalarGrid::AddInnerFace( int global_face_id, int bctype, int lcell, int rcell, int ftype )
{
	ReserveForAppend( this->global_faceid );
	ReserveForAppend( this->fBcTypes.data );
	ReserveForAppend( this->lc.data );
	ReserveForAppend( this->rc.data );
	ReserveForAppend( this->fTypes.data );

	this->global_faceid.push_back( global_face_id );
	this->fBcTypes.AddData( bctype );
	this->lc.AddData( lcell );
	this->rc.AddData( rcell );
	this->fTypes.AddData( ftype );
}

void ScalarGrid::AddInterface( int global_interface_id, int neighbor_zoneid, int neighbor_cellid )
{
	this->scalarIFace->AddInterface( global_interface_id, neighbor_zoneid, neighbor_cellid );
}

void ScalarGrid::ReconstructNode( const ScalarGrid & ggrid )
{
	const int nFaces = static_cast< int >( this->global_faceid.size() );
	const int nGlobalFaces = ggrid.GetNFaces();
	const int nGlobalNodes = ggrid.GetNNodes();
	if ( ggrid.xn.GetNElements() != static_cast< size_t >( nGlobalNodes ) ||
		 ggrid.yn.GetNElements() != static_cast< size_t >( nGlobalNodes ) ||
		 ggrid.zn.GetNElements() != static_cast< size_t >( nGlobalNodes ) ||
		 ggrid.faces.GetNElements() != static_cast< size_t >( nGlobalFaces ) )
	{
		throw std::runtime_error( "ScalarGrid::ReconstructNode: global node or face arrays have inconsistent sizes" );
	}

	std::vector< int > globalNodeIds;
	std::vector< std::vector< int > > reconstructedFaces;
	reconstructedFaces.reserve( nFaces );

	for ( int iFace = 0; iFace < nFaces; ++ iFace )
	{
		const int iGFace = this->global_faceid[ iFace ];
		if ( iGFace < 0 || iGFace >= nGlobalFaces )
		{
			throw std::runtime_error( "ScalarGrid::ReconstructNode: global face id is out of range" );
		}
		const std::vector< int > & face = ggrid.faces[ iGFace ];
		for ( const int globalNodeId : face )
		{
			if ( globalNodeId < 0 || globalNodeId >= nGlobalNodes )
			{
				throw std::runtime_error( "ScalarGrid::ReconstructNode: face references an invalid global node id" );
			}
			globalNodeIds.push_back( globalNodeId );
		}
		reconstructedFaces.push_back( face );
	}

	// Sort and deduplicate node ids to keep local numbering deterministic.
	std::sort( globalNodeIds.begin(), globalNodeIds.end() );
	globalNodeIds.erase( std::unique( globalNodeIds.begin(), globalNodeIds.end() ), globalNodeIds.end() );

	for ( std::vector< int > & face : reconstructedFaces )
	{
		for ( int & globalNodeId : face )
		{
			const auto localNode = std::lower_bound( globalNodeIds.begin(), globalNodeIds.end(), globalNodeId );
			globalNodeId = static_cast< int >( localNode - globalNodeIds.begin() );
		}
	}

	std::vector< Real > reconstructedX;
	std::vector< Real > reconstructedY;
	std::vector< Real > reconstructedZ;
	reconstructedX.reserve( globalNodeIds.size() );
	reconstructedY.reserve( globalNodeIds.size() );
	reconstructedZ.reserve( globalNodeIds.size() );
	for ( const int globalNodeId : globalNodeIds )
	{
		reconstructedX.push_back( ggrid.xn[ globalNodeId ] );
		reconstructedY.push_back( ggrid.yn[ globalNodeId ] );
		reconstructedZ.push_back( ggrid.zn[ globalNodeId ] );
	}

	// Commit reconstructed node and face data together, replacing any previous result.
	this->faces.data = std::move( reconstructedFaces );
	this->xn.data = std::move( reconstructedX );
	this->yn.data = std::move( reconstructedY );
	this->zn.data = std::move( reconstructedZ );

	// Face geometry, cell centers, and volumes are derived from the replaced mesh data.
	this->ResetGeometryData();
}
EndNameSpace

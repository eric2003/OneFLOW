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

#include "FaceSolver.h"
#include "ElementHome.h"
#include "Fatal.h"
#include "FaceTopo.h"
#include "BcRecord.h"
#include "CgnsSection.h"
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <unordered_set>


BeginNameSpace( ONEFLOW )

FaceSolver::FaceSolver()
{
    // [Refactored] faceBcNameIds, faceBcTypes, childFid are now value types, 
    // no need to 'new' them.
    this->faceTopo = std::make_unique< FaceTopo >();
}

FaceSolver::~FaceSolver()
    = default;

FaceTopo & FaceSolver::GetFaceTopo()
{
    return *this->faceTopo;
}

const FaceTopo & FaceSolver::GetFaceTopo() const
{
    return *this->faceTopo;
}

std::unique_ptr< FaceTopo > FaceSolver::ReleaseFaceTopo() noexcept
{
    return std::move( this->faceTopo );
}

bool FaceSolver::CheckBcFace( const IntSet & bcVertex, const IntField & nodeId ) const
{
    const int size = static_cast< int >( nodeId.size() );
    for ( int iNode = 0; iNode < size; ++ iNode )
    {
        if ( bcVertex.find( nodeId[ iNode ] ) == bcVertex.end() )
        {
            return false;
        }
    }
    return true;
}

void FaceSolver::StagePolygonFacesFromConnectivity(
    int nElement,
    const CgIntField & offsets,
    const CgIntField & connectivity,
    int faceType )
{
    if ( nElement < 0 ||
         static_cast< std::size_t >( nElement ) + 1 > offsets.size() )
    {
        throw std::runtime_error(
            "FaceSolver::StagePolygonFacesFromConnectivity: polygon offsets are incomplete" );
    }

    // Validate every span before mutating the face lookup or topology.
    for ( int iElem = 0; iElem < nElement; ++ iElem )
    {
        const CgInt start = offsets[ iElem ];
        const CgInt end = offsets[ iElem + 1 ];
        if ( start < 0 || end < start ||
             static_cast< std::size_t >( end ) > connectivity.size() )
        {
            throw std::runtime_error(
                "FaceSolver::StagePolygonFacesFromConnectivity: polygon connectivity span is invalid" );
        }
    }

    // Polygon faces are staged without an owning cell yet. Ownership is
    // assigned later via ClaimStagedFaceAsLeft. StageNewFace keeps topology
    // and BC slots aligned so StagedBcMatchesFaceCount() holds after scan.
    IntField faceNodes;
    for ( int iElem = 0; iElem < nElement; ++ iElem )
    {
        const CgInt start = offsets[ iElem ];
        const CgInt end = offsets[ iElem + 1 ];
        faceNodes.resize( 0 );
        for ( CgInt i = start; i < end; ++ i )
        {
            faceNodes.push_back( static_cast< int >( connectivity[ static_cast< std::size_t >( i ) ] ) );
        }

        auto [faceIndex, isNew] = faceLookup.FindOrAdd( faceNodes );
        if ( isNew )
        {
            this->StageNewFace( faceNodes, faceType, ONEFLOW::INVALID_INDEX );
        }
    }
}

void FaceSolver::ScanPolygonFace( CgnsSection & cgnsSection )
{
    this->StagePolygonFacesFromConnectivity(
        cgnsSection.nElement,
        cgnsSection.ePosList,
        cgnsSection.connList,
        cgnsSection.eType );
}

void FaceSolver::EnsureParallelFaceArraysSized()
{
    this->faceTopo->EnsureFaceSideArraysSized();
    const int nFaces = static_cast< int >( this->faceTopo->GetFaces().size() );
    this->faceBcTypes.resize( nFaces );
    this->faceBcNameIds.resize( nFaces );
    this->childFid.resize( nFaces );
}

void FaceSolver::ClaimPolyhedronCellsFromConnectivity(
    int nElement,
    const CgIntField & offsets,
    const CgIntField & connectivity )
{
    const auto & faces = this->faceTopo->GetFaces();

    if ( nElement < 0 ||
         static_cast< std::size_t >( nElement ) + 1 > offsets.size() )
    {
        throw std::runtime_error(
            "FaceSolver::ClaimPolyhedronCellsFromConnectivity: element offsets are incomplete" );
    }

    // Count existing cell incidences so a third cell cannot overwrite the right cell.
    const auto & faceFlags = this->faceTopo->GetFaceFlags();
    const auto & rightCells = this->faceTopo->GetRightCells();
    std::vector< int > faceUseCount( faces.size(), 0 );
    for ( std::size_t iFace = 0; iFace < faces.size(); ++ iFace )
    {
        if ( iFace < faceFlags.size() && faceFlags[ iFace ] != kFaceUnclaimed )
        {
            faceUseCount[ iFace ] =
                iFace < rightCells.size() && rightCells[ iFace ] != ONEFLOW::INVALID_INDEX ? 2 : 1;
        }
    }

    // Validate the complete section before changing cell adjacency.
    for ( int iElem = 0; iElem < nElement; ++ iElem )
    {
        const CgInt start = offsets[ iElem ];
        const CgInt end = offsets[ iElem + 1 ];
        if ( start < 0 || end < start ||
             static_cast< std::size_t >( end ) > connectivity.size() )
        {
            throw std::runtime_error(
                "FaceSolver::ClaimPolyhedronCellsFromConnectivity: connectivity span is invalid" );
        }

        std::unordered_set< std::size_t > referencedFaces;
        for ( CgInt i = start; i < end; ++ i )
        {
            const CgInt signedFaceId = connectivity[ static_cast< std::size_t >( i ) ];
            if ( signedFaceId == std::numeric_limits< CgInt >::min() )
            {
                throw std::runtime_error(
                    "FaceSolver::ClaimPolyhedronCellsFromConnectivity: face reference is out of range" );
            }

            const CgInt faceId = signedFaceId < 0 ? -signedFaceId : signedFaceId;
            if ( faceId < 0 || static_cast< std::size_t >( faceId ) >= faces.size() )
            {
                throw std::runtime_error(
                    "FaceSolver::ClaimPolyhedronCellsFromConnectivity: face reference is out of range" );
            }

            if ( ! referencedFaces.insert( static_cast< std::size_t >( faceId ) ).second )
            {
                throw std::runtime_error(
                    "FaceSolver::ClaimPolyhedronCellsFromConnectivity: polyhedron references the same face more than once" );
            }

            const std::size_t polygonFaceId = static_cast< std::size_t >( faceId );
            if ( ++ faceUseCount[ polygonFaceId ] > 2 )
            {
                throw std::runtime_error(
                    "FaceSolver::ClaimPolyhedronCellsFromConnectivity: face is referenced by more than two cells" );
            }
        }
    }

    for ( int iElem = 0; iElem < nElement; ++ iElem )
    {
        const CgInt start = offsets[ iElem ];
        const CgInt end = offsets[ iElem + 1 ];

        for ( CgInt i = start; i < end; ++ i )
        {
            const CgInt signedFaceId = connectivity[ static_cast< std::size_t >( i ) ];
            const CgInt faceId = signedFaceId < 0 ? -signedFaceId : signedFaceId;
            const int polygonFaceId = static_cast< int >( faceId );
            const int faceFlag =
                this->faceTopo->GetFaceFlags()[ static_cast< std::size_t >( polygonFaceId ) ];

            if ( faceFlag == kFaceUnclaimed )
            {
                this->ClaimStagedFaceAsLeft( polygonFaceId, iElem );
            }
            else
            {
                this->AttachStagedFaceAsRight( polygonFaceId, iElem );
            }
        }
    }
}

void FaceSolver::ScanPolyhedronElement( CgnsSection & cgnsSection )
{
    this->ClaimPolyhedronCellsFromConnectivity(
        cgnsSection.nElement,
        cgnsSection.ePosList,
        cgnsSection.connList );
}


void FaceSolver::StageNewFace( const IntField & faceNodes, int faceType, int leftCellId )
{
    // Single write path for a newly discovered unique face. All parallel
    // arrays (topology, flags, cell sides, face-local positions, BC slots)
    // grow together so callers never need ad-hoc push_back or EnsureFaceSideArraysSized.
    this->faceTopo->GetFaces().push_back( faceNodes );
    this->faceTopo->GetFaceTypes().push_back( faceType );
    this->faceTopo->GetFaceFlags().push_back( kFaceUnclaimed );
    this->faceTopo->GetLeftCells().push_back( leftCellId );
    this->faceTopo->GetRightCells().push_back( ONEFLOW::INVALID_INDEX );
    this->faceTopo->lPosition.push_back( 0 );
    this->faceTopo->rPosition.push_back( 0 );

    this->faceBcTypes.push_back( ONEFLOW::INVALID_INDEX );
    this->faceBcNameIds.push_back( ONEFLOW::INVALID_INDEX );
    this->childFid.push_back( IntField{} );
}

void FaceSolver::SetStagedFaceBc( int faceId, int bcType, int bcNameId )
{
    if ( faceId < 0 ||
         static_cast< size_t >( faceId ) >= this->faceBcTypes.size() ||
         static_cast< size_t >( faceId ) >= this->faceBcNameIds.size() )
    {
        throw std::out_of_range( "FaceSolver::SetStagedFaceBc: face ID is out of range" );
    }
    this->faceBcTypes[ faceId ] = bcType;
    this->faceBcNameIds[ faceId ] = bcNameId;
}

void FaceSolver::ClaimStagedFaceAsLeft( int faceId, int cellId )
{
    if ( faceId < 0 ||
         static_cast< size_t >( faceId ) >= this->faceTopo->GetFaces().size() ||
         static_cast< size_t >( faceId ) >= this->faceTopo->GetFaceFlags().size() )
    {
        throw std::out_of_range( "FaceSolver::ClaimStagedFaceAsLeft: face ID is out of range" );
    }
    if ( this->faceTopo->GetFaceFlags()[ static_cast< size_t >( faceId ) ] != kFaceUnclaimed )
    {
        throw std::logic_error( "FaceSolver::ClaimStagedFaceAsLeft: face is already claimed" );
    }
    if ( ! this->StagedBcMatchesFaceCount() )
    {
        throw std::logic_error(
            "FaceSolver::ClaimStagedFaceAsLeft: staged face arrays are not aligned" );
    }

    // StageNewFace already sized topology, positions, and BC slots. Claim only
    // updates ownership; it must not resize.
    this->faceTopo->GetFaceFlags()[ static_cast< size_t >( faceId ) ] = kFaceClaimed;
    this->faceTopo->GetLeftCells()[ static_cast< size_t >( faceId ) ] = cellId;
    this->faceTopo->GetRightCells()[ static_cast< size_t >( faceId ) ] = ONEFLOW::INVALID_INDEX;
    this->SetStagedFaceBc( faceId, ONEFLOW::INVALID_INDEX, ONEFLOW::INVALID_INDEX );
}

void FaceSolver::AttachStagedFaceAsRight( int faceId, int cellId )
{
    if ( faceId < 0 ||
         static_cast< size_t >( faceId ) >= this->faceTopo->GetFaces().size() ||
         static_cast< size_t >( faceId ) >= this->faceTopo->GetFaceFlags().size() )
    {
        throw std::out_of_range( "FaceSolver::AttachStagedFaceAsRight: face ID is out of range" );
    }
    if ( this->faceTopo->GetFaceFlags()[ static_cast< size_t >( faceId ) ] == kFaceUnclaimed )
    {
        throw std::logic_error( "FaceSolver::AttachStagedFaceAsRight: face is not claimed yet" );
    }

    this->faceTopo->GetRightCells()[ static_cast< size_t >( faceId ) ] = cellId;
}

void FaceSolver::BindCellToExistingFace( int faceId, int cellId )
{
    if ( faceId < 0 ||
         static_cast< size_t >( faceId ) >= this->faceTopo->GetFaces().size() ||
         static_cast< size_t >( faceId ) >= this->faceTopo->GetLeftCells().size() ||
         static_cast< size_t >( faceId ) >= this->faceTopo->GetRightCells().size() )
    {
        throw std::out_of_range( "FaceSolver::BindCellToExistingFace: face ID is out of range" );
    }

    auto & leftCells = this->faceTopo->GetLeftCells();
    auto & rightCells = this->faceTopo->GetRightCells();
    const std::size_t id = static_cast< std::size_t >( faceId );

    if ( leftCells[ id ] == ONEFLOW::INVALID_INDEX )
    {
        // Face was staged without an owner (or left cleared during split).
        leftCells[ id ] = cellId;
        return;
    }
    if ( rightCells[ id ] == ONEFLOW::INVALID_INDEX && leftCells[ id ] != cellId )
    {
        rightCells[ id ] = cellId;
    }
}

bool FaceSolver::StagedBcMatchesFaceCount() const
{
    // Full parallel-array contract after StageNewFace: topology, ownership
    // slots, face-local positions, BC attributes, and child-face lists.
    const std::size_t nFaces = this->faceTopo->GetFaces().size();
    return this->faceBcTypes.size() == nFaces &&
           this->faceBcNameIds.size() == nFaces &&
           this->childFid.size() == nFaces &&
           this->faceTopo->GetFaceTypes().size() == nFaces &&
           this->faceTopo->GetFaceFlags().size() == nFaces &&
           this->faceTopo->GetLeftCells().size() == nFaces &&
           this->faceTopo->GetRightCells().size() == nFaces &&
           this->faceTopo->lPosition.size() == nFaces &&
           this->faceTopo->rPosition.size() == nFaces;
}

void FaceSolver::ScanElementFace( CgIntField & eNodeId, int eType, int eId )
{
    const UnitElement & elementDefinition = ElementHome::GetUnitElementDefinition( eType );

    // Composite elements are not involved in this analysis.
    const int elementFaceCount = elementDefinition.GetElementFaceCount();
    for ( int localFaceIndex = 0; localFaceIndex < elementFaceCount; ++ localFaceIndex )
    {
        const IntField & localFaceNodeIndices =
            elementDefinition.GetElementFace( localFaceIndex );
        const int faceType = elementDefinition.GetFaceType( localFaceIndex );
         
        int nNodes = localFaceNodeIndices.size();

        IntField aNodeId;
        for ( int iNode = 0; iNode < nNodes; ++ iNode )
        {
            aNodeId.push_back( eNodeId[ localFaceNodeIndices[ iNode ] ] );
        }                                                              

        auto [gFid, isNew] = this->faceLookup.FindOrAdd( aNodeId );

        if ( isNew )
        {
            this->StageNewFace( aNodeId, faceType, eId );
        }
        else
        {
            this->BindCellToExistingFace( gFid, eId );
        }
    }
}

int FaceSolver::AssignBcToMatchingFaces(
    const IntSet & bcVertex,
    int bcType,
    int bcNameId,
    bool boundaryFacesOnly )
{
    const int nFaces = static_cast< int >( this->faceTopo->GetLeftCells().size() );
    int nAssigned = 0;

    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        if ( this->faceBcTypes[ iFace ] != ONEFLOW::INVALID_INDEX )
        {
            continue;
        }
        if ( boundaryFacesOnly &&
             this->faceTopo->GetRightCells()[ iFace ] != ONEFLOW::INVALID_INDEX )
        {
            continue;
        }
        if ( this->CheckBcFace( bcVertex, this->faceTopo->GetFaces()[ iFace ] ) )
        {
            this->SetStagedFaceBc( iFace, bcType, bcNameId );
            ++ nAssigned;
        }
    }
    return nAssigned;
}

void FaceSolver::ScanBcFace( IntSet & bcVertex, int bcType, int bcNameId )
{
    // Boundary faces only: open right cell and no BC yet.
    this->AssignBcToMatchingFaces( bcVertex, bcType, bcNameId, /*boundaryFacesOnly=*/true );
}

void FaceSolver::ScanBcFaceDetail( IntSet & bcVertex, int bcType, int bcNameId )
{
    // Match by vertex set only; do not require an open right cell.
    this->AssignBcToMatchingFaces( bcVertex, bcType, bcNameId, /*boundaryFacesOnly=*/false );
}

void FaceSolver::ScanInterfaceBc()
{
    const int nFaces = static_cast< int >( this->faceTopo->GetLeftCells().size() );
    const int bcNameId = -1;

    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        if ( this->faceBcTypes[ iFace ] == ONEFLOW::INVALID_INDEX )
        {
            this->SetStagedFaceBc( iFace, BCTypeNull, bcNameId );
        }
    }
}


int FaceSolver::AppendBoundaryBcToRecord( BcRecord & bcRecord ) const
{
    if ( ! this->StagedBcMatchesFaceCount() )
    {
        throw std::logic_error(
            "FaceSolver::AppendBoundaryBcToRecord: staged face arrays are not aligned" );
    }

    const auto & rightCells = this->faceTopo->GetRightCells();
    const int nFaces = static_cast< int >( rightCells.size() );
    int nBFaces = 0;

    for ( int iFace = 0; iFace < nFaces; ++ iFace )
    {
        if ( rightCells[ iFace ] == ONEFLOW::INVALID_INDEX )
        {
            bcRecord.bcType.push_back( this->faceBcTypes[ iFace ] );
            bcRecord.bcNameId.push_back( this->faceBcNameIds[ iFace ] );
            ++ nBFaces;
        }
    }
    return nBFaces;
}

int FaceSolver::GetNSimpleFace()
{
    int nSimpleFace = 0;

    for ( int iFace = 0; iFace < this->faceTopo->GetFaces().size(); ++ iFace )
    {
        int nCFace = this->childFid[ iFace ].size();
        if ( nCFace == 0 )
        {
            ++ nSimpleFace;
        }
    }
    return nSimpleFace;
}


EndNameSpace

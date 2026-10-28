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

#include "FaceSearch.h"
#include "PointLocator.h"
#include "Grid.h"
#include "NodeMesh.h"
#include "IFaceLink.h"
#include "Dimension.h"
#include "Fatal.h"
#include <algorithm>
#include <stdexcept>


BeginNameSpace( ONEFLOW )

namespace
{

void ResolveRegisteredChildFaces(
    const FaceSearchBasic & faceSearch,
    const LinkField & candidateFaces,
    IntField & childFaceIds,
    std::vector< int > & candidateIndices )
{
    IntField resolvedChildFaceIds;
    std::vector< int > resolvedCandidateIndices;
    for ( int iCandidate = 0;
         iCandidate < static_cast< int >( candidateFaces.size() );
         ++ iCandidate )
    {
        const int childFaceId = faceSearch.FindFace( candidateFaces[ iCandidate ] );
        if ( childFaceId == ONEFLOW::INVALID_INDEX ) continue;

        if ( std::find(
                 resolvedChildFaceIds.begin(), resolvedChildFaceIds.end(), childFaceId ) !=
             resolvedChildFaceIds.end() )
        {
            throw std::logic_error(
                "ResolveRegisteredChildFaces: duplicate child face for parent" );
        }

        resolvedChildFaceIds.push_back( childFaceId );
        resolvedCandidateIndices.push_back( iCandidate );
    }

    childFaceIds.swap( resolvedChildFaceIds );
    candidateIndices.swap( resolvedCandidateIndices );
}

}

int FaceSearchBasic::AddFace( const IntField & faceNode )
{
    // HXLookup automatically sorts the nodes for uniqueness
    auto [faceId, isNew] = lookup_.FindOrAdd( faceNode );

    if ( isNew )
    {
        // New face: store the original node order
        faceNodeRings.push_back( faceNode );
    }

    return faceId;   // always return a valid ID
}

int FaceSearchBasic::FindFace( const IntField & faceNode ) const
{
    return lookup_.Find( faceNode );   // returns INVALID_INDEX if not found
}

void FaceSearchBasic::Clear()
{
    faceNodeRings.clear();
    lookup_.Clear();
}

const IntField & FaceSearchBasic::GetFaceNodeRing( int faceId ) const
{
    if ( faceId < 0 || static_cast< size_t >( faceId ) >= this->faceNodeRings.size() )
    {
        throw std::out_of_range( "FaceSearchBasic::GetFaceNodeRing: face ID is out of range" );
    }
    return this->faceNodeRings[ static_cast< size_t >( faceId ) ];
}


void InterfaceFaceSplitTable::ClearParentMapping( int parentFaceId )
{
    if ( parentFaceId < 0 ||
         static_cast< std::size_t >( parentFaceId ) >= this->childFaceIdsOfParent.size() )
    {
        throw std::out_of_range(
            "InterfaceFaceSplitTable::ClearParentMapping: parent face ID is out of range" );
    }

    const bool hasParentSpecificTables =
        ! this->childFaceNodeIdsByParent.empty() ||
        ! this->childFaceNodeFlagsByParent.empty();
    if ( hasParentSpecificTables &&
         ( static_cast< std::size_t >( parentFaceId ) >= this->childFaceNodeIdsByParent.size() ||
           static_cast< std::size_t >( parentFaceId ) >= this->childFaceNodeFlagsByParent.size() ) )
    {
        throw std::logic_error(
            "InterfaceFaceSplitTable::ClearParentMapping: parent-specific tables are misaligned" );
    }

    this->childFaceIdsOfParent[ parentFaceId ].clear();
    if ( hasParentSpecificTables )
    {
        this->childFaceNodeIdsByParent[ parentFaceId ].clear();
        this->childFaceNodeFlagsByParent[ parentFaceId ].clear();
    }
}

void InterfaceFaceSplitTable::SetChildMapping(
    int parentFaceId,
    int childFaceId,
    const IntField & localNodeIds,
    const IntField & localNodeFlags )
{
    if ( parentFaceId < 0 ||
         static_cast< std::size_t >( parentFaceId ) >= this->childFaceIdsOfParent.size() ||
         childFaceId < 0 ||
         static_cast< std::size_t >( childFaceId ) >= this->childFaceNodeIds.size() ||
         static_cast< std::size_t >( childFaceId ) >= this->childFaceNodeFlags.size() )
    {
        throw std::out_of_range(
            "InterfaceFaceSplitTable::SetChildMapping: face ID is out of range" );
    }
    if ( localNodeIds.size() != localNodeFlags.size() )
    {
        throw std::logic_error(
            "InterfaceFaceSplitTable::SetChildMapping: node IDs and flags have different sizes" );
    }

    if ( static_cast< std::size_t >( parentFaceId ) >= this->childFaceNodeIdsByParent.size() ||
         static_cast< std::size_t >( parentFaceId ) >= this->childFaceNodeFlagsByParent.size() )
    {
        throw std::logic_error(
            "InterfaceFaceSplitTable::SetChildMapping: parent-specific tables are misaligned" );
    }

    const IntField & parentChildren = this->childFaceIdsOfParent[ parentFaceId ];
    const auto childPosition = std::find(
        parentChildren.begin(), parentChildren.end(), childFaceId );
    if ( childPosition == parentChildren.end() )
    {
        throw std::logic_error(
            "InterfaceFaceSplitTable::SetChildMapping: child is not registered for parent" );
    }
    if ( std::find( childPosition + 1, parentChildren.end(), childFaceId ) !=
         parentChildren.end() )
    {
        throw std::logic_error(
            "InterfaceFaceSplitTable::SetChildMapping: child is registered more than once for parent" );
    }

    const std::size_t parentIndex = static_cast< std::size_t >( parentFaceId );
    const std::size_t childIndex = static_cast< std::size_t >(
        childPosition - parentChildren.begin() );
    LinkField & parentNodeIds = this->childFaceNodeIdsByParent[ parentIndex ];
    LinkField & parentNodeFlags = this->childFaceNodeFlagsByParent[ parentIndex ];
    if ( parentNodeIds.size() != parentNodeFlags.size() ||
         parentNodeIds.size() > parentChildren.size() )
    {
        throw std::logic_error(
            "InterfaceFaceSplitTable::SetChildMapping: parent mapping rows are misaligned" );
    }

    if ( childIndex < parentNodeIds.size() &&
         ( ! parentNodeIds[ childIndex ].empty() ||
           ! parentNodeFlags[ childIndex ].empty() ) &&
         ( parentNodeIds[ childIndex ] != localNodeIds ||
           parentNodeFlags[ childIndex ] != localNodeFlags ) )
    {
        throw std::logic_error(
            "InterfaceFaceSplitTable::SetChildMapping: conflicting mapping for parent child" );
    }

    parentNodeIds.resize( parentChildren.size() );
    parentNodeFlags.resize( parentChildren.size() );
    parentNodeIds[ childIndex ] = localNodeIds;
    parentNodeFlags[ childIndex ] = localNodeFlags;

    if ( this->childFaceNodeIds[ childFaceId ].empty() &&
         this->childFaceNodeFlags[ childFaceId ].empty() )
    {
        this->childFaceNodeIds[ childFaceId ] = localNodeIds;
        this->childFaceNodeFlags[ childFaceId ] = localNodeFlags;
    }
}

void FaceSearch::Clear()
{
    FaceSearchBasic::Clear();
    this->interfaceSplit.Clear();
}

void FaceSearch::CalcNewFaceId( IFaceLink & iFaceLink )
{
    const int nFaces = static_cast< int >( this->Size() );
    this->interfaceSplit.ResetForFaceCount( static_cast< std::size_t >( nFaces ) );

    if ( Dim::dimension == ONEFLOW::THREE_D )
    {
        for ( int iFace = 0; iFace < nFaces; ++ iFace )
        {
            this->SplitQuad2Tri( iFace );
        }
    }
    else if ( Dim::dimension == ONEFLOW::TWO_D )
    {
        for ( int iFace = 0; iFace < nFaces; ++ iFace )
        {
            this->SplitLine( iFaceLink, iFace );
        }
    }
}

void FaceSearch::SplitQuad2Tri( int faceId )
{
    const IntField & nodeId = this->GetFaceNodeRing( faceId );

    // A direct call also invalidates stale results when this face is not a quad.
    this->interfaceSplit.ClearParentMapping( faceId );

    if ( nodeId.size() != 4 ) return;

    LinkField localTriId, localTriFlag;
    BuildQuadDiagonalLocalTris( localTriId, localTriFlag );

    LinkField triId;
    MapLocalTrisToGlobalNodes( nodeId, localTriId, triId );

    // Resolve child identities before publishing parent-relative mappings.
    IntField childFaceIds;
    std::vector< int > childCandidateIndices;
    ResolveRegisteredChildFaces( *this, triId, childFaceIds, childCandidateIndices );

    IntField & registeredChildFaceIds = this->interfaceSplit.childFaceIdsOfParent[ faceId ];
    for ( std::size_t iChild = 0; iChild < childFaceIds.size(); ++ iChild )
    {
        const int childFaceId = childFaceIds[ iChild ];
        const int candidateIndex = childCandidateIndices[ iChild ];
        registeredChildFaceIds.push_back( childFaceId );
        this->interfaceSplit.SetChildMapping(
            faceId, childFaceId, localTriId[ candidateIndex ], localTriFlag[ candidateIndex ] );
    }
}

void FaceSearch::SplitLine( IFaceLink & iFaceLink, int globalFaceId )
{
    const IntField & nodeId = this->GetFaceNodeRing( globalFaceId );

    // Rebuild one parent's mapping instead of accumulating stale child IDs.
    this->interfaceSplit.ClearParentMapping( globalFaceId );

    if ( nodeId.size() != 2 ) return;

    LinkField localLineId, localLineFlag;
    LinkField lineId;
    if ( ! this->TrySplitLineAtMidpoint(
             iFaceLink, globalFaceId, nodeId, localLineId, localLineFlag, lineId ) )
    {
        return;
    }

    // Resolve child identities before publishing parent-relative mappings.
    IntField childFaceIds;
    std::vector< int > childCandidateIndices;
    ResolveRegisteredChildFaces( *this, lineId, childFaceIds, childCandidateIndices );

    IntField & registeredChildFaceIds =
        this->interfaceSplit.childFaceIdsOfParent[ globalFaceId ];
    for ( std::size_t iChild = 0; iChild < childFaceIds.size(); ++ iChild )
    {
        const int childFaceId = childFaceIds[ iChild ];
        const int candidateIndex = childCandidateIndices[ iChild ];
        registeredChildFaceIds.push_back( childFaceId );
        this->interfaceSplit.SetChildMapping(
            globalFaceId, childFaceId,
            localLineId[ candidateIndex ], localLineFlag[ candidateIndex ] );
    }
}

void BuildLineSplitMapping(
    int firstNode,
    int secondNode,
    int midpointLocalNode,
    int midpointGlobalNode,
    LinkField & localLineId,
    LinkField & localLineFlag,
    LinkField & lineId )
{
    localLineId.clear();
    localLineFlag.clear();
    lineId.clear();

    localLineId.push_back( IntField{ 0, midpointLocalNode } );
    localLineId.push_back( IntField{ midpointLocalNode, 1 } );

    localLineFlag.push_back( IntField{ 1, -1 } );
    localLineFlag.push_back( IntField{ -1, 1 } );

    lineId.push_back( IntField{ firstNode, midpointGlobalNode } );
    lineId.push_back( IntField{ midpointGlobalNode, secondNode } );
}

void BuildQuadDiagonalLocalTris( LinkField & localTriId, LinkField & localTriFlag )
{
    localTriId.clear();
    localTriFlag.clear();

    // Diagonal 0-2
    localTriId.push_back( IntField{ 0, 1, 2 } );
    localTriId.push_back( IntField{ 2, 3, 0 } );
    localTriFlag.push_back( IntField{ 1, 1, 1 } );
    localTriFlag.push_back( IntField{ 1, 1, 1 } );

    // Diagonal 1-3
    localTriId.push_back( IntField{ 3, 0, 1 } );
    localTriId.push_back( IntField{ 1, 2, 3 } );
    localTriFlag.push_back( IntField{ 1, 1, 1 } );
    localTriFlag.push_back( IntField{ 1, 1, 1 } );
}

void MapLocalTrisToGlobalNodes(
    const IntField & nodeId,
    const LinkField & localTriId,
    LinkField & triId )
{
    // Validate the complete candidate topology before replacing the output.
    for ( std::size_t iTri = 0; iTri < localTriId.size(); ++ iTri )
    {
        if ( localTriId[ iTri ].size() != 3 )
        {
            throw std::logic_error(
                "MapLocalTrisToGlobalNodes: each local triangle must have three corners" );
        }

        for ( int localCorner : localTriId[ iTri ] )
        {
            if ( localCorner < 0 ||
                 static_cast< std::size_t >( localCorner ) >= nodeId.size() )
            {
                throw std::out_of_range(
                    "MapLocalTrisToGlobalNodes: local corner index is out of range" );
            }
        }
    }

    LinkField mappedTriIds( localTriId.size() );
    for ( std::size_t iTri = 0; iTri < localTriId.size(); ++ iTri )
    {
        for ( int localCorner : localTriId[ iTri ] )
        {
            mappedTriIds[ iTri ].push_back( nodeId[ localCorner ] );
        }
    }

    triId.swap( mappedTriIds );
}

bool FaceSearch::TrySplitLineAtMidpoint(
    IFaceLink & iFaceLink,
    int globalFaceId,
    const IntField & nodeId,
    LinkField & localLineId,
    LinkField & localLineFlag,
    LinkField & lineId )
{
    if ( ! iFaceLink.point_search )
    {
        throw std::logic_error(
            "FaceSearch::TrySplitLineAtMidpoint: point locator is not initialized" );
    }
    PointLocator * point_search = iFaceLink.point_search.get();
    const int p0 = nodeId[ 0 ];
    const int p1 = nodeId[ 1 ];
    RealField coor0( 3 ), coor1( 3 ), coor2( 3 );
    point_search->GetPoint( p0, coor0[ 0 ], coor0[ 1 ], coor0[ 2 ] );
    point_search->GetPoint( p1, coor1[ 0 ], coor1[ 1 ], coor1[ 2 ] );
    for ( int m = 0; m < 3; ++ m )
    {
        coor2[ m ] = half * ( coor0[ m ] + coor1[ m ] );
    }
    const int midpointGlobalNode =
        point_search->FindPoint( coor2[ 0 ], coor2[ 1 ], coor2[ 2 ] );
    if ( midpointGlobalNode == ONEFLOW::INVALID_INDEX ) return false;

    const LinkField & globalFaceZones = iFaceLink.gI2Zid;
    if ( globalFaceId < 0 ||
         static_cast< size_t >( globalFaceId ) >= globalFaceZones.size() )
    {
        throw std::out_of_range(
            "FaceSearch::TrySplitLineAtMidpoint: global face ID is out of range" );
    }

    const IntField & zoneIds = globalFaceZones[ globalFaceId ];
    if ( zoneIds.size() != 1 )
    {
        throw std::logic_error(
            "FaceSearch::TrySplitLineAtMidpoint: a split line must reference exactly one zone" );
    }

    const int zoneIndex = zoneIds[ 0 ];
    Grid & grid = iFaceLink.GetGrid( zoneIndex );
    int pId = grid.nodeMesh->FindPoint(
        coor2[ 0 ], coor2[ 1 ], coor2[ 2 ], point_search->GetTol() );
    if ( pId == ONEFLOW::INVALID_INDEX )
    {
        pId = static_cast< int >( grid.nodeMesh->GetNumberOfNodes() );
        grid.nodeMesh->AddPoint( coor2[ 0 ], coor2[ 1 ], coor2[ 2 ] );
    }
    // The midpoint was found above; reuse its global ID without another lookup.
    BuildLineSplitMapping(
        p0, p1, pId, midpointGlobalNode, localLineId, localLineFlag, lineId );

    return true;
}

EndNameSpace

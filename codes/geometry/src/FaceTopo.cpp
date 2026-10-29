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

#include "FaceTopo.h"
#include "BcRecord.h"
#include "Boundary.h"
#include "IFaceLink.h"
#include "InterFace.h"
#include "FaceSearch.h"
#include "Grid.h"
#include "NodeMesh.h"
#include "HXMath.h"
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>
#include <utility>


BeginNameSpace( ONEFLOW )

namespace
{
int GetBoundaryFaceId(
    const Grid & grid, int interfaceFaceId, const IntField & cellIds,
    const char * operation )
{
    if ( ! grid.interFace )
    {
        throw std::logic_error( std::string( operation ) + ": interface data is not initialized" );
    }
    const IntField & interfaceToBoundary = grid.interFace->i2b;
    if ( interfaceFaceId < 0 ||
         static_cast< size_t >( interfaceFaceId ) >= interfaceToBoundary.size() )
    {
        throw std::out_of_range( std::string( operation ) + ": interface face ID is out of range" );
    }

    const int boundaryFaceId = interfaceToBoundary[ interfaceFaceId ];
    if ( boundaryFaceId < 0 ||
         static_cast< size_t >( boundaryFaceId ) >= cellIds.size() )
    {
        throw std::out_of_range( std::string( operation ) + ": boundary face ID is out of range" );
    }
    return boundaryFaceId;
}

} // namespace

const IntField & GetInterfaceChildFaceIds(
    const LinkField & zoneToLocalToGlobal,
    const InterfaceFaceSplitTable & splitTable,
    int zoneId,
    int localInterfaceFaceId )
{
    if ( zoneId < 0 ||
         static_cast< size_t >( zoneId ) >= zoneToLocalToGlobal.size() ||
         localInterfaceFaceId < 0 ||
         static_cast< size_t >( localInterfaceFaceId ) >=
             zoneToLocalToGlobal[ zoneId ].size() )
    {
        throw std::out_of_range(
            "GetInterfaceChildFaceIds: local interface face is out of range" );
    }

    const int globalFaceId = zoneToLocalToGlobal[ zoneId ][ localInterfaceFaceId ];
    const LinkField & childFaces = splitTable.childFaceIdsOfParent;
    if ( globalFaceId < 0 ||
         static_cast< size_t >( globalFaceId ) >= childFaces.size() )
    {
        throw std::out_of_range(
            "GetInterfaceChildFaceIds: global interface face is out of range" );
    }

    const IntField & childFaceIds = childFaces[ globalFaceId ];

    // This accessor resolves parent-to-child identity only. The split table
    // does not carry the global face count, so child ID bounds are checked by
    // consumers that have access to the global face registry.
    return childFaceIds;
}

const IntField & GetInterfaceChildFaceIds(
    const IFaceLink & interfaceLink, int zoneId, int localInterfaceFaceId )
{
    return GetInterfaceChildFaceIds(
        interfaceLink.GetLocalToGlobalMap(),
        interfaceLink.GetInterfaceSplit(),
        zoneId,
        localInterfaceFaceId );
}

FaceTopo::FaceTopo()
{
    this->bcManager = std::make_unique< BcManager >();
}

void FaceTopo::BindGrid( Grid & grid )
{
    this->grid = &grid;
}

Grid & FaceTopo::GetGrid()
{
    if ( ! this->grid )
    {
        throw std::logic_error( "FaceTopo::GetGrid: grid is not bound" );
    }
    return *this->grid;
}

const Grid & FaceTopo::GetGrid() const
{
    if ( ! this->grid )
    {
        throw std::logic_error( "FaceTopo::GetGrid: grid is not bound" );
    }
    return *this->grid;
}

BcRecord & FaceTopo::GetBcRecord()
{
    return *this->bcManager->bcRecord;
}

const BcRecord & FaceTopo::GetBcRecord() const
{
    return *this->bcManager->bcRecord;
}

IntField & FaceTopo::GetFaceFlags()
{
    return this->faceFlags;
}

const IntField & FaceTopo::GetFaceFlags() const
{
    return this->faceFlags;
}

IntField & FaceTopo::GetLeftCells()
{
    return this->lCells;
}

const IntField & FaceTopo::GetLeftCells() const
{
    return this->lCells;
}

IntField & FaceTopo::GetRightCells()
{
    return this->rCells;
}

const IntField & FaceTopo::GetRightCells() const
{
    return this->rCells;
}

LinkField & FaceTopo::GetFaces()
{
    return this->faces;
}

const LinkField & FaceTopo::GetFaces() const
{
    return this->faces;
}

IntField & FaceTopo::GetFaceTypes()
{
    return this->fTypes;
}

const IntField & FaceTopo::GetFaceTypes() const
{
    return this->fTypes;
}

void FaceTopo::PrepareBoundaryConditions()
{
    this->bcManager->PreProcess();
}

bool FaceTopo::HasInterfaceBoundary() const
{
    return this->bcManager->ExistInterface();
}

FaceTopo::~FaceTopo() = default;

HXSize_t FaceTopo::CalcTotalFaceNodes() const
{
    HXSize_t totalNumFaceNodes = 0;
    HXSize_t nFaces = this->GetFaceCount();
    for ( HXSize_t iFace = 0; iFace < nFaces; ++ iFace )
    {
        totalNumFaceNodes += faces[ iFace ].size();
    }
    return totalNumFaceNodes;
}

HXSize_t FaceTopo::GetBoundaryFaceCount() const
{
    return this->bcManager->bcRecord->GetBoundaryFaceCount();
}

HXSize_t FaceTopo::GetNBFaces() const
{
    return GetBoundaryFaceCount();
}

void FaceTopo::SetBoundaryFaceCount( HXSize_t boundaryFaceCount )
{
    this->bcManager->bcRecord->Init( boundaryFaceCount );
}

void FaceTopo::SetNBFaces( HXSize_t nBFaces )
{
    SetBoundaryFaceCount( nBFaces );
}

void FaceTopo::EnsureFaceSideArraysSized()
{
    const std::size_t nFaces = this->faces.size();
    this->lCells.resize( nFaces );
    this->rCells.resize( nFaces );
    this->lPosition.resize( nFaces );
    this->rPosition.resize( nFaces );
    // Keep claim flags aligned with faces when present or when growing from empty.
    this->faceFlags.resize( nFaces );
}

void FaceTopo::ModifyFaceNodeId( IFaceLink & iFaceLink )
{
    this->SetNewFace2Node( iFaceLink );
    this->SetNewFace2Cell( iFaceLink );
}


IntField ResolveChildFaceNodeIds(
    const IntField & parentFaceNodes,
    const IntField & childNodeIds,
    const IntField & childNodeFlags,
    bool hasNodeMesh,
    HXSize_t nodeMeshCount )
{
    if ( childNodeIds.size() != childNodeFlags.size() )
    {
        throw std::logic_error(
            "ResolveChildFaceNodeIds: child node IDs and flags have different sizes" );
    }

    IntField resolved;
    resolved.reserve( childNodeIds.size() );

    for ( size_t iNode = 0; iNode < childNodeIds.size(); ++ iNode )
    {
        const int flag = childNodeFlags[ iNode ];
        const int mappedNodeId = childNodeIds[ iNode ];
        int nodeIndex = 0;

        if ( flag == 1 )
        {
            if ( mappedNodeId < 0 ||
                 static_cast< size_t >( mappedNodeId ) >= parentFaceNodes.size() )
            {
                throw std::out_of_range(
                    "ResolveChildFaceNodeIds: local child node ID is out of range" );
            }
            nodeIndex = parentFaceNodes[ mappedNodeId ];
        }
        else
        {
            if ( ! hasNodeMesh ||
                 mappedNodeId < 0 ||
                 static_cast< size_t >( mappedNodeId ) >= nodeMeshCount )
            {
                throw std::out_of_range(
                    "ResolveChildFaceNodeIds: absolute child node ID is out of range" );
            }
            nodeIndex = mappedNodeId;
        }
        resolved.push_back( nodeIndex );
    }
    return resolved;
}

void FaceTopo::SetNewFace2Node( IFaceLink & iFaceLink )
{
    const HXSize_t nFaces = this->GetFaceCount();
    const HXSize_t nBoundaryFaces = this->bcManager->bcRecord->GetBoundaryFaceCount();

    // A face has one node list and one face type; boundary faces must be a subset.
    if ( this->faces.size() != nFaces || nBoundaryFaces > nFaces )
    {
        throw std::logic_error(
            "FaceTopo::SetNewFace2Node: face, face-type, and boundary-face counts are inconsistent" );
    }

    const int nBFaces = static_cast< int >( nBoundaryFaces );
    const int nFacesInt = static_cast< int >( nFaces );
    const IntField & bcType = this->bcManager->bcRecord->bcType;

    LinkField childGlobalIdsPerInterfaceLocal;
    CollectBoundaryInterfaceChildFaceIds(
        bcType,
        nBFaces,
        iFaceLink.GetLocalToGlobalMap(),
        iFaceLink.MutableInterfaceSplit(),
        this->GetGrid().id,
        childGlobalIdsPerInterfaceLocal );

    const bool hasNodeMesh = static_cast< bool >( this->GetGrid().nodeMesh );
    const HXSize_t nodeMeshCount = hasNodeMesh
        ? this->GetGrid().nodeMesh->GetNumberOfNodes()
        : 0;

    std::vector< LinkField > childNodeIdsPerInterfaceLocal;
    std::vector< LinkField > childNodeFlagsPerInterfaceLocal;
    CollectBoundaryInterfaceChildNodeMappings(
        bcType,
        nBFaces,
        iFaceLink.GetLocalToGlobalMap(),
        iFaceLink.MutableInterfaceSplit(),
        this->GetGrid().id,
        childNodeIdsPerInterfaceLocal,
        childNodeFlagsPerInterfaceLocal );

    RebuildFaceNodeListsAfterInterfaceSplit(
        this->faces,
        bcType,
        nBFaces,
        nFacesInt,
        childGlobalIdsPerInterfaceLocal,
        childNodeIdsPerInterfaceLocal,
        childNodeFlagsPerInterfaceLocal,
        hasNodeMesh,
        nodeMeshCount,
        this->rebuildFaces );
}

void FaceTopo::SetNewFace2Cell( IFaceLink & iFaceLink )
{
    const int nBFaces = this->bcManager->bcRecord->GetBoundaryFaceCount();
    const int nFaces = static_cast< int >( this->GetFaceCount() );
    const IntField & bcType = this->bcManager->bcRecord->bcType;

    IntField childCountPerInterfaceLocal;
    CollectBoundaryInterfaceChildCounts(
        bcType,
        nBFaces,
        iFaceLink.GetLocalToGlobalMap(),
        iFaceLink.MutableInterfaceSplit(),
        this->GetGrid().id,
        childCountPerInterfaceLocal );

    RebuildFaceCellSidesAfterInterfaceSplit(
        bcType,
        this->lCells,
        this->rCells,
        nBFaces,
        nFaces,
        childCountPerInterfaceLocal,
        this->GetGrid().nCells,
        this->rebuildLeftCells,
        this->rebuildRightCells );
}


void FaceTopo::ModifyBoundaryInformation( IFaceLink & iFaceLink )
{
    const int zoneId = this->GetGrid().id;

    // Live parent globals and split bookkeeping go through IFaceLink so zone
    // bounds and nChild sizing stay next to the dual-map tables.
    const IntField & parentGlobalFaceIds = iFaceLink.GetL2gForZone( zoneId );
    const int nIFaces = static_cast< int >( parentGlobalFaceIds.size() );

    LinkField childGlobalIdsPerParent;
    childGlobalIdsPerParent.resize( static_cast< size_t >( nIFaces ) );
    for ( int iFid = 0; iFid < nIFaces; ++ iFid )
    {
        childGlobalIdsPerParent[ iFid ] = GetInterfaceChildFaceIds(
            iFaceLink.GetLocalToGlobalMap(),
            iFaceLink.MutableInterfaceSplit(),
            zoneId,
            iFid );
    }

    iFaceLink.EnsureSplitBookkeepingMatchesZones();

    RebuildZoneL2gNewAfterInterfaceSplit(
        parentGlobalFaceIds,
        childGlobalIdsPerParent,
        iFaceLink.MutableL2gNewForZone( zoneId ),
        iFaceLink.MutableNChildForZone( zoneId ) );

    this->ResetNumberOfBoundaryCondition( iFaceLink );
}

void FaceTopo::ResetNumberOfBoundaryCondition( IFaceLink & iFaceLink )
{
    const int nBFaces = this->bcManager->bcRecord->GetBoundaryFaceCount();
    const IntField & bcType = this->bcManager->bcRecord->bcType;
    const IntField & bcNameId = this->bcManager->bcRecord->bcNameId;

    // Child counts in boundary-interface order (same walk as ExpandBcRecord...).
    IntField childCountPerInterfaceLocal;
    CollectBoundaryInterfaceChildCounts(
        bcType,
        nBFaces,
        iFaceLink.GetLocalToGlobalMap(),
        iFaceLink.MutableInterfaceSplit(),
        this->GetGrid().id,
        childCountPerInterfaceLocal );

    ExpandBcRecordForSplitInterfaces(
        bcType,
        bcNameId,
        childCountPerInterfaceLocal,
        this->bcManager->bcRecordNew->bcType,
        this->bcManager->bcRecordNew->bcNameId );

    this->ConstructNewInterfaceMap( iFaceLink );
}

void FaceTopo::ConstructNewInterfaceMap( IFaceLink & iFaceLink )
{
    const int zoneId = this->GetGrid().id;

    BuildInverseInterfaceMapForZone(
        zoneId,
        iFaceLink.MutableL2gNewForZone( zoneId ),
        iFaceLink.MutableGlobalToZoneIdsNew(),
        iFaceLink.MutableGlobalToLocalIdsNew() );
}

void FaceTopo::UpdateOtherTopologyTerm()
{
    this->bcManager->Update();
    this->lCells = this->rebuildLeftCells;
    this->rCells = this->rebuildRightCells;
    this->faces = this->rebuildFaces;

    // Drop staging storage so a completed rebuild cannot leak stale buffers.
    this->ClearRebuildBuffers();

    // The face IDs may have changed when interface faces were split.
    this->c2f.clear();
}

void FaceTopo::ClearRebuildBuffers()
{
    this->rebuildFaces.clear();
    this->rebuildLeftCells.clear();
    this->rebuildRightCells.clear();
}

bool FaceTopo::HasPendingTopologyRebuild() const
{
    return ! this->rebuildFaces.empty() ||
           ! this->rebuildLeftCells.empty() ||
           ! this->rebuildRightCells.empty();
}

void FaceTopo::GenerateI2B( InterFace & interFace )
{
    this->bcManager->bcRecord->GenerateI2B( interFace );
}

bool FaceTopo::GetSId( int iFace, int iPosition, int & sId )
{
    const Grid & grid = this->GetGrid();
    const int boundaryFaceId = GetBoundaryFaceId(
        grid, iFace, this->lCells, "FaceTopo::GetSId" );
    sId = this->lCells[ boundaryFaceId ];

    return true;
}

bool FaceTopo::GetTId( int iFace, int iPosition, int & tId )
{
    const Grid & grid = this->GetGrid();
    const int boundaryFaceId = GetBoundaryFaceId(
        grid, iFace, this->rCells, "FaceTopo::GetTId" );
    tId = this->rCells[ boundaryFaceId ];

    return true;
}

void FaceTopo::CalcC2C( LinkField & c2c )
{
    const Grid & grid = this->GetGrid();
    const int nCells = grid.nCells;
    const HXSize_t nBFaces = this->GetBoundaryFaceCount();
    const HXSize_t nFaces = this->GetFaceCount();
    const IntField & leftCells = this->GetLeftCells();
    const IntField & rightCells = this->GetRightCells();
    const IntField & boundaryTypes = this->bcManager->bcRecord->bcType;

    if ( nCells < 0 || nBFaces > nFaces ||
         leftCells.size() != nFaces || rightCells.size() != nFaces ||
         boundaryTypes.size() != nBFaces )
    {
        throw std::logic_error(
            "FaceTopo::CalcC2C: inconsistent face-cell topology sizes" );
    }

    LinkField nextC2c( static_cast< HXSize_t >( nCells ) );

    // Interface boundaries may reference ghost cells outside the local cell range.
    for ( HXSize_t iFace = 0; iFace < nBFaces; ++ iFace )
    {
        const int leftCell = leftCells[ iFace ];
        if ( leftCell < 0 || leftCell >= nCells )
        {
            throw std::out_of_range(
                "FaceTopo::CalcC2C: boundary face left cell ID is out of range" );
        }

        if ( BC::IsInterfaceBc( boundaryTypes[ iFace ] ) )
        {
            nextC2c[ leftCell ].push_back( rightCells[ iFace ] );
        }
    }

    for ( HXSize_t iFace = nBFaces; iFace < nFaces; ++ iFace )
    {
        const int leftCell = leftCells[ iFace ];
        const int rightCell = rightCells[ iFace ];
        if ( leftCell < 0 || leftCell >= nCells ||
             rightCell < 0 || rightCell >= nCells )
        {
            throw std::out_of_range(
                "FaceTopo::CalcC2C: internal face cell ID is out of range" );
        }
        if ( leftCell == rightCell )
        {
            throw std::logic_error(
                "FaceTopo::CalcC2C: internal face must connect distinct cells" );
        }

        nextC2c[ leftCell ].push_back( rightCell );
        nextC2c[ rightCell ].push_back( leftCell );
    }

    c2c = nextC2c;
}

IntField BuildBoundaryFirstPermutation( const IntField & rightCells )
{
    const HXSize_t nFaces = rightCells.size();
    IntField newToOld( nFaces );
    HXSize_t nextFace = 0;

    for ( HXSize_t iFace = 0; iFace < nFaces; ++ iFace )
    {
        if ( rightCells[ iFace ] == INVALID_INDEX )
        {
            newToOld[ nextFace ++ ] = static_cast< int >( iFace );
        }
    }
    for ( HXSize_t iFace = 0; iFace < nFaces; ++ iFace )
    {
        if ( rightCells[ iFace ] != INVALID_INDEX )
        {
            newToOld[ nextFace ++ ] = static_cast< int >( iFace );
        }
    }
    return newToOld;
}


HXSize_t CountBoundaryFaces( const IntField & rightCells )
{
    HXSize_t nBoundary = 0;
    for ( const int rightCell : rightCells )
    {
        if ( rightCell == INVALID_INDEX )
        {
            ++ nBoundary;
        }
    }
    return nBoundary;
}

void ValidateBoundaryBcLayout( const FaceTopo & faceTopo )
{
    const IntField & rightCells = faceTopo.GetRightCells();
    const HXSize_t nFaces = rightCells.size();
    const HXSize_t nBoundary = CountBoundaryFaces( rightCells );
    const HXSize_t nBcRecord = static_cast< HXSize_t >(
        faceTopo.GetBcRecord().GetBoundaryFaceCount() );

    if ( nBcRecord != nBoundary )
    {
        throw std::logic_error(
            "ValidateBoundaryBcLayout: BcRecord size does not match boundary face count" );
    }

    for ( HXSize_t iFace = 0; iFace < nBoundary; ++ iFace )
    {
        if ( rightCells[ iFace ] != INVALID_INDEX )
        {
            throw std::logic_error(
                "ValidateBoundaryBcLayout: face in boundary prefix has an interior right cell" );
        }
    }
    for ( HXSize_t iFace = nBoundary; iFace < nFaces; ++ iFace )
    {
        if ( rightCells[ iFace ] == INVALID_INDEX )
        {
            throw std::logic_error(
                "ValidateBoundaryBcLayout: boundary face found after boundary prefix" );
        }
    }
}

void ReorderBoundaryFirstAndValidate( FaceTopo & faceTopo )
{
    faceTopo.ReorderLink();
    ValidateBoundaryBcLayout( faceTopo );
}

void FaceTopo::ReorderLink()
{
    const HXSize_t nFaces = this->GetFaceCount();
    const IntField & rightCells = this->GetRightCells();

    if ( this->faces.size() != nFaces ||
         this->GetLeftCells().size() != nFaces ||
         rightCells.size() != nFaces )
    {
        throw std::logic_error( "FaceTopo::ReorderLink: inconsistent face topology sizes" );
    }

    const auto validateOptionalFaceData = [ nFaces ]( const IntField & values, const char * name )
    {
        if ( ! values.empty() && values.size() != nFaces )
        {
            throw std::logic_error(
                std::string( "FaceTopo::ReorderLink: inconsistent " ) + name + " size" );
        }
    };
    validateOptionalFaceData( this->lPosition, "left face-position" );
    validateOptionalFaceData( this->rPosition, "right face-position" );
    validateOptionalFaceData( this->faceFlags, "face-flag" );

    const IntField newToOld = BuildBoundaryFirstPermutation( rightCells );

    const auto reorderFaceData = [ &newToOld, nFaces ]( auto & values )
    {
        if ( values.empty() )
        {
            return;
        }
        const auto original = values;
        for ( HXSize_t iFace = 0; iFace < nFaces; ++ iFace )
        {
            values[ iFace ] = original[ newToOld[ iFace ] ];
        }
    };

    reorderFaceData( this->faces );
    reorderFaceData( this->fTypes );
    reorderFaceData( this->lCells );
    reorderFaceData( this->rCells );
    reorderFaceData( this->lPosition );
    reorderFaceData( this->rPosition );
    reorderFaceData( this->faceFlags );

    // Cell-to-face adjacency stores face indices and must be rebuilt after reordering.
    this->c2f.clear();
}


void ExpandBcRecordForSplitInterfaces(
    const IntField & bcType,
    const IntField & bcNameId,
    const IntField & childCountPerInterfaceLocal,
    IntField & outBcType,
    IntField & outBcNameId )
{
    if ( bcType.size() != bcNameId.size() )
    {
        throw std::logic_error(
            "ExpandBcRecordForSplitInterfaces: bcType and bcNameId size mismatch" );
    }

    outBcType.clear();
    outBcNameId.clear();

    int interfaceLocal = 0;
    for ( std::size_t iFace = 0; iFace < bcType.size(); ++ iFace )
    {
        const int bc = bcType[ iFace ];
        const int nameId = bcNameId[ iFace ];
        if ( BC::IsInterfaceBc( bc ) )
        {
            if ( interfaceLocal < 0 ||
                 static_cast< std::size_t >( interfaceLocal ) >=
                     childCountPerInterfaceLocal.size() )
            {
                throw std::out_of_range(
                    "ExpandBcRecordForSplitInterfaces: interface child-count index out of range" );
            }
            const int nCFace = childCountPerInterfaceLocal[ interfaceLocal ];
            ++ interfaceLocal;
            const int nEmit = ( nCFace > 0 ) ? nCFace : 1;
            for ( int i = 0; i < nEmit; ++ i )
            {
                outBcType.push_back( bc );
                outBcNameId.push_back( nameId );
            }
        }
        else
        {
            outBcType.push_back( bc );
            outBcNameId.push_back( nameId );
        }
    }
}


void BuildInverseInterfaceMapForZone(
    int zoneId,
    const IntField & l2gNewForZone,
    LinkField & gI2ZidNew,
    LinkField & g2lNew )
{
    if ( zoneId < 0 )
    {
        throw std::out_of_range(
            "BuildInverseInterfaceMapForZone: zone ID is negative" );
    }
    if ( gI2ZidNew.size() != g2lNew.size() )
    {
        throw std::logic_error(
            "BuildInverseInterfaceMapForZone: gI2ZidNew and g2lNew size mismatch" );
    }

    // Validate all input IDs and existing rows before changing either table.
    // One zone must not register the same global face more than once: otherwise
    // the duplicate could be mistaken for a remote participant during lookup.
    std::unordered_set< int > seenGlobalFaceIds;
    seenGlobalFaceIds.reserve( l2gNewForZone.size() );
    for ( const int gIFid : l2gNewForZone )
    {
        if ( gIFid < 0 )
        {
            throw std::out_of_range(
                "BuildInverseInterfaceMapForZone: global face id is negative" );
        }
        if ( ! seenGlobalFaceIds.insert( gIFid ).second )
        {
            throw std::logic_error(
                "BuildInverseInterfaceMapForZone: duplicate global face ID in zone mapping" );
        }

        if ( static_cast< std::size_t >( gIFid ) < gI2ZidNew.size() &&
             gI2ZidNew[ gIFid ].size() != g2lNew[ gIFid ].size() )
        {
            throw std::logic_error(
                "BuildInverseInterfaceMapForZone: per-face zone/local lists differ in size" );
        }
    }

    for ( int localIFid = 0; localIFid < static_cast< int >( l2gNewForZone.size() ); ++ localIFid )
    {
        const int gIFid = l2gNewForZone[ localIFid ];
        const std::size_t need = static_cast< std::size_t >( gIFid ) + 1;
        if ( gI2ZidNew.size() < need )
        {
            gI2ZidNew.resize( need );
            g2lNew.resize( need );
        }

        // Rebuild buffers start with the live inverse maps. Refresh this
        // zone's slot when the global face already exists; otherwise a stale
        // local index or duplicate self-participant would survive the rebuild.
        IntField & zoneIds = gI2ZidNew[ gIFid ];
        IntField & localFaceIds = g2lNew[ gIFid ];
        bool foundZone = false;
        for ( std::size_t slot = 0; slot < zoneIds.size(); ++ slot )
        {
            if ( zoneIds[ slot ] == zoneId )
            {
                localFaceIds[ slot ] = localIFid;
                foundZone = true;
                break;
            }
        }

        // Keep the legacy two-participant capacity. A third distinct zone is
        // ignored, but a zone already in the row must be refreshed above.
        if ( ! foundZone && zoneIds.size() < 2 )
        {
            zoneIds.push_back( zoneId );
            localFaceIds.push_back( localIFid );
        }
    }
}


void RebuildZoneL2gNewAfterInterfaceSplit(
    const IntField & parentGlobalFaceIds,
    const LinkField & childGlobalIdsPerParent,
    IntField & l2gNewForZone,
    IntField & nChildForZone )
{
    if ( parentGlobalFaceIds.size() != childGlobalIdsPerParent.size() )
    {
        throw std::logic_error(
            "RebuildZoneL2gNewAfterInterfaceSplit: parent and child tables differ in size" );
    }

    // Validate the complete replacement map before changing either output.
    // Negative global IDs cannot identify a valid face and must not leak into
    // the candidate mapping for later inverse-map construction.
    for ( std::size_t iFid = 0; iFid < parentGlobalFaceIds.size(); ++ iFid )
    {
        if ( parentGlobalFaceIds[ iFid ] < 0 )
        {
            throw std::out_of_range(
                "RebuildZoneL2gNewAfterInterfaceSplit: parent global face ID is negative" );
        }
        for ( const int childGlobalFaceId : childGlobalIdsPerParent[ iFid ] )
        {
            if ( childGlobalFaceId < 0 )
            {
                throw std::out_of_range(
                    "RebuildZoneL2gNewAfterInterfaceSplit: child global face ID is negative" );
            }
        }
    }

    l2gNewForZone.clear();
    nChildForZone.assign( parentGlobalFaceIds.size(), 0 );

    for ( std::size_t iFid = 0; iFid < parentGlobalFaceIds.size(); ++ iFid )
    {
        const IntField & childFaceIds = childGlobalIdsPerParent[ iFid ];
        const int gFid = parentGlobalFaceIds[ iFid ];
        const int nCFace = static_cast< int >( childFaceIds.size() );

        if ( nCFace > 0 )
        {
            nChildForZone[ iFid ] = nCFace;
            for ( int iCFace = 0; iCFace < nCFace; ++ iCFace )
            {
                l2gNewForZone.push_back( childFaceIds[ iCFace ] );
                nChildForZone.push_back( 0 );
            }
        }
        else
        {
            l2gNewForZone.push_back( gFid );
        }
    }
}


void RebuildFaceCellSidesAfterInterfaceSplit(
    const IntField & bcType,
    const IntField & leftCells,
    const IntField & rightCells,
    int nBFaces,
    int nFaces,
    const IntField & childCountPerInterfaceLocal,
    int nCells,
    IntField & outLeftCells,
    IntField & outRightCells )
{
    if ( nBFaces < 0 || nFaces < 0 || nBFaces > nFaces )
    {
        throw std::logic_error(
            "RebuildFaceCellSidesAfterInterfaceSplit: invalid face counts" );
    }
    if ( static_cast< std::size_t >( nBFaces ) > bcType.size() ||
         static_cast< std::size_t >( nFaces ) > leftCells.size() ||
         static_cast< std::size_t >( nFaces ) > rightCells.size() )
    {
        throw std::logic_error(
            "RebuildFaceCellSidesAfterInterfaceSplit: arrays shorter than face counts" );
    }

    outLeftCells.clear();
    outRightCells.clear();

    int interfaceLocal = 0;
    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int nEmit = 1;
        if ( BC::IsInterfaceBc( bcType[ iFace ] ) )
        {
            if ( interfaceLocal < 0 ||
                 static_cast< std::size_t >( interfaceLocal ) >=
                     childCountPerInterfaceLocal.size() )
            {
                throw std::out_of_range(
                    "RebuildFaceCellSidesAfterInterfaceSplit: interface child-count index out of range" );
            }
            const int nCFace = childCountPerInterfaceLocal[ interfaceLocal ];
            ++ interfaceLocal;
            nEmit = ( nCFace > 0 ) ? nCFace : 1;
        }
        for ( int i = 0; i < nEmit; ++ i )
        {
            outLeftCells.push_back( leftCells[ iFace ] );
        }
    }

    const int nBFaceNew = static_cast< int >( outLeftCells.size() );
    for ( int iFace = 0; iFace < nBFaceNew; ++ iFace )
    {
        outRightCells.push_back( iFace + nCells );
    }

    for ( int iFace = nBFaces; iFace < nFaces; ++ iFace )
    {
        outLeftCells.push_back( leftCells[ iFace ] );
        outRightCells.push_back( rightCells[ iFace ] );
    }
}


void CollectBoundaryInterfaceChildFaceIds(
    const IntField & bcType,
    int nBFaces,
    const LinkField & zoneToLocalToGlobal,
    const InterfaceFaceSplitTable & splitTable,
    int zoneId,
    LinkField & outChildGlobalIdsPerInterfaceLocal )
{
    if ( nBFaces < 0 || static_cast< std::size_t >( nBFaces ) > bcType.size() )
    {
        throw std::logic_error(
            "CollectBoundaryInterfaceChildFaceIds: invalid boundary face count" );
    }

    outChildGlobalIdsPerInterfaceLocal.clear();
    int interfaceLocal = 0;
    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        if ( ! BC::IsInterfaceBc( bcType[ iFace ] ) )
        {
            continue;
        }
        outChildGlobalIdsPerInterfaceLocal.push_back(
            GetInterfaceChildFaceIds(
                zoneToLocalToGlobal, splitTable, zoneId, interfaceLocal ) );
        ++ interfaceLocal;
    }
}

void CollectBoundaryInterfaceChildNodeMappings(
    const IntField & bcType,
    int nBFaces,
    const LinkField & zoneToLocalToGlobal,
    const InterfaceFaceSplitTable & splitTable,
    int zoneId,
    std::vector< LinkField > & outNodeIdsPerInterfaceLocal,
    std::vector< LinkField > & outNodeFlagsPerInterfaceLocal )
{
    if ( nBFaces < 0 || static_cast< std::size_t >( nBFaces ) > bcType.size() )
    {
        throw std::logic_error(
            "CollectBoundaryInterfaceChildNodeMappings: invalid boundary face count" );
    }
    if ( zoneId < 0 || static_cast< std::size_t >( zoneId ) >= zoneToLocalToGlobal.size() )
    {
        throw std::out_of_range(
            "CollectBoundaryInterfaceChildNodeMappings: zone ID is out of range" );
    }

    // Build both outputs off to the side. A malformed later parent must not
    // leave callers with a partial set of node mappings.
    std::vector< LinkField > nodeIdsByInterfaceLocal;
    std::vector< LinkField > nodeFlagsByInterfaceLocal;
    const bool hasParentSpecificTables =
        ! splitTable.childFaceNodeIdsByParent.empty() ||
        ! splitTable.childFaceNodeFlagsByParent.empty();
    bool anyParentSpecificMapping = false;
    if ( hasParentSpecificTables )
    {
        if ( splitTable.childFaceNodeIdsByParent.size() !=
             splitTable.childFaceNodeFlagsByParent.size() )
        {
            throw std::logic_error(
                "CollectBoundaryInterfaceChildNodeMappings: parent-specific tables are misaligned" );
        }
        for ( std::size_t parentIndex = 0;
              parentIndex < splitTable.childFaceNodeIdsByParent.size();
              ++ parentIndex )
        {
            if ( ! splitTable.childFaceNodeIdsByParent[ parentIndex ].empty() ||
                 ! splitTable.childFaceNodeFlagsByParent[ parentIndex ].empty() )
            {
                anyParentSpecificMapping = true;
                break;
            }
        }
    }

    int interfaceLocal = 0;
    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        if ( ! BC::IsInterfaceBc( bcType[ iFace ] ) ) continue;

        const IntField & childIds = GetInterfaceChildFaceIds(
            zoneToLocalToGlobal, splitTable, zoneId, interfaceLocal );
        const int globalParentId = zoneToLocalToGlobal[ zoneId ][ interfaceLocal ];
        if ( childIds.empty() )
        {
            // An unsplit interface face needs no node mapping. Keep the empty
            // row so the rebuild stage can preserve its parent face unchanged.
            nodeIdsByInterfaceLocal.emplace_back();
            nodeFlagsByInterfaceLocal.emplace_back();
            ++ interfaceLocal;
            continue;
        }

        if ( hasParentSpecificTables )
        {
            if ( globalParentId < 0 ||
                 static_cast< std::size_t >( globalParentId ) >= splitTable.childFaceNodeIdsByParent.size() )
            {
                throw std::logic_error(
                    "CollectBoundaryInterfaceChildNodeMappings: parent-specific tables do not contain the parent" );
            }

            const LinkField & nodeIds = splitTable.childFaceNodeIdsByParent[ globalParentId ];
            const LinkField & nodeFlags = splitTable.childFaceNodeFlagsByParent[ globalParentId ];
            if ( ! nodeIds.empty() || ! nodeFlags.empty() )
            {
                if ( nodeIds.size() != childIds.size() ||
                     nodeFlags.size() != childIds.size() )
                {
                    throw std::logic_error(
                        "CollectBoundaryInterfaceChildNodeMappings: parent-specific rows do not match child count" );
                }
                nodeIdsByInterfaceLocal.push_back( nodeIds );
                nodeFlagsByInterfaceLocal.push_back( nodeFlags );
                ++ interfaceLocal;
                continue;
            }

            // Once any parent-specific mapping has been published, an empty
            // row for a split parent is missing data, not a legacy fallback.
            if ( anyParentSpecificMapping )
            {
                throw std::logic_error(
                    "CollectBoundaryInterfaceChildNodeMappings: missing parent-specific mapping for split parent" );
            }
        }

        // Older callers and serialized test fixtures may leave all
        // parent-specific rows empty and populate only child-global rows.
        // Use the compatibility representation only when no parent-specific
        // mapping has been published anywhere.
        LinkField nodeIdsByChild;
        LinkField nodeFlagsByChild;
        for ( const int childFaceId : childIds )
        {
            if ( childFaceId < 0 ||
                 static_cast< std::size_t >( childFaceId ) >= splitTable.childFaceNodeIds.size() ||
                 static_cast< std::size_t >( childFaceId ) >= splitTable.childFaceNodeFlags.size() )
            {
                throw std::out_of_range(
                    "CollectBoundaryInterfaceChildNodeMappings: child face ID is out of range" );
            }
            const IntField & nodeIds = splitTable.childFaceNodeIds[ childFaceId ];
            const IntField & nodeFlags = splitTable.childFaceNodeFlags[ childFaceId ];
            if ( nodeIds.size() != nodeFlags.size() )
            {
                throw std::logic_error(
                    "CollectBoundaryInterfaceChildNodeMappings: child node IDs and flags have different sizes" );
            }
            nodeIdsByChild.push_back( nodeIds );
            nodeFlagsByChild.push_back( nodeFlags );
        }
        nodeIdsByInterfaceLocal.push_back( nodeIdsByChild );
        nodeFlagsByInterfaceLocal.push_back( nodeFlagsByChild );
        ++ interfaceLocal;
    }

    outNodeIdsPerInterfaceLocal.swap( nodeIdsByInterfaceLocal );
    outNodeFlagsPerInterfaceLocal.swap( nodeFlagsByInterfaceLocal );
}

void CollectBoundaryInterfaceChildCounts(
    const IntField & bcType,
    int nBFaces,
    const LinkField & zoneToLocalToGlobal,
    const InterfaceFaceSplitTable & splitTable,
    int zoneId,
    IntField & outChildCountPerInterfaceLocal )
{
    LinkField childGlobalIdsPerInterfaceLocal;
    CollectBoundaryInterfaceChildFaceIds(
        bcType,
        nBFaces,
        zoneToLocalToGlobal,
        splitTable,
        zoneId,
        childGlobalIdsPerInterfaceLocal );

    outChildCountPerInterfaceLocal.clear();
    outChildCountPerInterfaceLocal.reserve( childGlobalIdsPerInterfaceLocal.size() );
    for ( std::size_t i = 0; i < childGlobalIdsPerInterfaceLocal.size(); ++ i )
    {
        outChildCountPerInterfaceLocal.push_back(
            static_cast< int >( childGlobalIdsPerInterfaceLocal[ i ].size() ) );
    }
}


void RebuildFaceNodeListsAfterInterfaceSplit(
    const LinkField & faces,
    const IntField & bcType,
    int nBFaces,
    int nFaces,
    const LinkField & childGlobalIdsPerInterfaceLocal,
    const LinkField & childFaceNodeIds,
    const LinkField & childFaceNodeFlags,
    bool hasNodeMesh,
    HXSize_t nodeMeshCount,
    LinkField & outFaces )
{
    if ( nBFaces < 0 || nFaces < 0 || nBFaces > nFaces )
    {
        throw std::logic_error(
            "RebuildFaceNodeListsAfterInterfaceSplit: invalid face counts" );
    }
    if ( static_cast< std::size_t >( nFaces ) > faces.size() ||
         static_cast< std::size_t >( nBFaces ) > bcType.size() )
    {
        throw std::logic_error(
            "RebuildFaceNodeListsAfterInterfaceSplit: arrays shorter than face counts" );
    }

    // Adapt the legacy child-global tables to the parent-local row layout.
    // The second overload owns the single face-list reconstruction algorithm.
    std::vector< LinkField > childNodeIdsPerInterfaceLocal;
    std::vector< LinkField > childNodeFlagsPerInterfaceLocal;
    int interfaceLocal = 0;
    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        if ( ! BC::IsInterfaceBc( bcType[ iFace ] ) ) continue;
        if ( static_cast< std::size_t >( interfaceLocal ) >=
             childGlobalIdsPerInterfaceLocal.size() )
        {
            throw std::out_of_range(
                "RebuildFaceNodeListsAfterInterfaceSplit: interface child list index out of range" );
        }

        const IntField & childFaceIds =
            childGlobalIdsPerInterfaceLocal[ interfaceLocal ];
        ++ interfaceLocal;
        LinkField nodeIdsByChild;
        LinkField nodeFlagsByChild;
        for ( const int childFaceId : childFaceIds )
        {
            if ( childFaceId < 0 ||
                 static_cast< std::size_t >( childFaceId ) >= childFaceNodeIds.size() ||
                 static_cast< std::size_t >( childFaceId ) >= childFaceNodeFlags.size() )
            {
                throw std::out_of_range(
                    "RebuildFaceNodeListsAfterInterfaceSplit: child face id out of range" );
            }
            nodeIdsByChild.push_back( childFaceNodeIds[ childFaceId ] );
            nodeFlagsByChild.push_back( childFaceNodeFlags[ childFaceId ] );
        }
        childNodeIdsPerInterfaceLocal.push_back( nodeIdsByChild );
        childNodeFlagsPerInterfaceLocal.push_back( nodeFlagsByChild );
    }

    RebuildFaceNodeListsAfterInterfaceSplit(
        faces,
        bcType,
        nBFaces,
        nFaces,
        childGlobalIdsPerInterfaceLocal,
        childNodeIdsPerInterfaceLocal,
        childNodeFlagsPerInterfaceLocal,
        hasNodeMesh,
        nodeMeshCount,
        outFaces );
}

void RebuildFaceNodeListsAfterInterfaceSplit(
    const LinkField & faces,
    const IntField & bcType,
    int nBFaces,
    int nFaces,
    const LinkField & childGlobalIdsPerInterfaceLocal,
    const std::vector< LinkField > & childNodeIdsPerInterfaceLocal,
    const std::vector< LinkField > & childNodeFlagsPerInterfaceLocal,
    bool hasNodeMesh,
    HXSize_t nodeMeshCount,
    LinkField & outFaces )
{
    if ( nBFaces < 0 || nFaces < 0 || nBFaces > nFaces )
    {
        throw std::logic_error(
            "RebuildFaceNodeListsAfterInterfaceSplit: invalid face counts" );
    }
    if ( static_cast< std::size_t >( nFaces ) > faces.size() ||
         static_cast< std::size_t >( nBFaces ) > bcType.size() )
    {
        throw std::logic_error(
            "RebuildFaceNodeListsAfterInterfaceSplit: arrays shorter than face counts" );
    }
    if ( childNodeIdsPerInterfaceLocal.size() != childGlobalIdsPerInterfaceLocal.size() ||
         childNodeFlagsPerInterfaceLocal.size() != childGlobalIdsPerInterfaceLocal.size() )
    {
        throw std::logic_error(
            "RebuildFaceNodeListsAfterInterfaceSplit: parent mapping count differs from child list count" );
    }

    LinkField rebuiltFaces;
    int interfaceLocal = 0;
    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        if ( ! BC::IsInterfaceBc( bcType[ iFace ] ) )
        {
            rebuiltFaces.push_back( faces[ iFace ] );
            continue;
        }
        if ( static_cast< std::size_t >( interfaceLocal ) >= childGlobalIdsPerInterfaceLocal.size() )
        {
            throw std::out_of_range(
                "RebuildFaceNodeListsAfterInterfaceSplit: interface child list index out of range" );
        }

        const IntField & childFaceIds = childGlobalIdsPerInterfaceLocal[ interfaceLocal ];
        const LinkField & childNodeIds = childNodeIdsPerInterfaceLocal[ interfaceLocal ];
        const LinkField & childNodeFlags = childNodeFlagsPerInterfaceLocal[ interfaceLocal ];
        ++ interfaceLocal;

        if ( childFaceIds.empty() )
        {
            rebuiltFaces.push_back( faces[ iFace ] );
            continue;
        }
        if ( childNodeIds.size() != childFaceIds.size() ||
             childNodeFlags.size() != childFaceIds.size() )
        {
            throw std::logic_error(
                "RebuildFaceNodeListsAfterInterfaceSplit: parent mapping rows do not match child IDs" );
        }
        for ( std::size_t childIndex = 0; childIndex < childFaceIds.size(); ++ childIndex )
        {
            if ( childNodeIds[ childIndex ].size() != childNodeFlags[ childIndex ].size() )
            {
                throw std::logic_error(
                    "RebuildFaceNodeListsAfterInterfaceSplit: child node IDs and flags have different sizes" );
            }
            rebuiltFaces.push_back(
                ResolveChildFaceNodeIds(
                    faces[ iFace ],
                    childNodeIds[ childIndex ],
                    childNodeFlags[ childIndex ],
                    hasNodeMesh,
                    nodeMeshCount ) );
        }
    }
    for ( int iFace = nBFaces; iFace < nFaces; ++ iFace )
    {
        rebuiltFaces.push_back( faces[ iFace ] );
    }
    outFaces.swap( rebuiltFaces );
}



void RunSymmetricInterfaceSplitMappingPipeline(
    const LinkField & l2gIn,
    const LinkField & gI2ZidIn,
    const LinkField & g2lIn,
    const IntField & parentGlobalFaceIds,
    const LinkField & childGlobalIdsPerParent,
    const IntField & bcType,
    const IntField & bcNameId,
    const IntField & childCountPerInterfaceLocal,
    InterfaceSplitMappingPipelineResult & out )
{
    if ( l2gIn.empty() )
    {
        throw std::invalid_argument(
            "RunSymmetricInterfaceSplitMappingPipeline: l2gIn must not be empty" );
    }
    if ( childGlobalIdsPerParent.size() != parentGlobalFaceIds.size() )
    {
        throw std::invalid_argument(
            "RunSymmetricInterfaceSplitMappingPipeline: parent/child table size mismatch" );
    }

    // Build a complete candidate result. A failed stage must not publish a
    // mixture of new child counts / BC records and old interface maps.
    InterfaceSplitMappingPipelineResult next;
    LinkField gI2ZidNew, g2lNew, l2gNew;
    InitNewInterfaceMappingBuffers(
        gI2ZidIn, g2lIn, l2gIn, gI2ZidNew, g2lNew, l2gNew );

    next.nChild.resize( l2gIn.size() );

    for ( std::size_t zid = 0; zid < l2gIn.size(); ++ zid )
    {
        RebuildZoneL2gNewAfterInterfaceSplit(
            parentGlobalFaceIds,
            childGlobalIdsPerParent,
            l2gNew[ zid ],
            next.nChild[ zid ] );
    }

    for ( std::size_t zid = 0; zid < l2gIn.size(); ++ zid )
    {
        BuildInverseInterfaceMapForZone(
            static_cast< int >( zid ),
            l2gNew[ zid ],
            gI2ZidNew,
            g2lNew );
    }

    ExpandBcRecordForSplitInterfaces(
        bcType, bcNameId, childCountPerInterfaceLocal,
        next.outBcType, next.outBcNameId );

    next.l2g = l2gIn;
    next.gI2Zid = gI2ZidIn;
    next.g2l = g2lIn;
    CommitInterfaceMappingBuffers(
        next.gI2Zid, next.g2l, next.l2g,
        gI2ZidNew, g2lNew, l2gNew,
        "RunSymmetricInterfaceSplitMappingPipeline" );

    out = std::move( next );
}


void CollectUnreferencedGlobalInterfaceFaces(
    const LinkField & localToGlobal,
    std::size_t globalFaceRowCount,
    IntField & outUnreferencedSorted )
{
    // Validate fully before mutating the output (strong exception guarantee).
    for ( const IntField & row : localToGlobal )
    {
        for ( const int gid : row )
        {
            if ( gid < 0 )
            {
                throw std::out_of_range(
                    "CollectUnreferencedGlobalInterfaceFaces: negative global face id" );
            }
            if ( static_cast< std::size_t >( gid ) >= globalFaceRowCount )
            {
                throw std::out_of_range(
                    "CollectUnreferencedGlobalInterfaceFaces: global face id exceeds row count" );
            }
        }
    }

    outUnreferencedSorted.clear();
    if ( globalFaceRowCount == 0 )
    {
        return;
    }

    std::vector< char > referenced( globalFaceRowCount, 0 );
    for ( const IntField & row : localToGlobal )
    {
        for ( const int gid : row )
        {
            referenced[ static_cast< std::size_t >( gid ) ] = 1;
        }
    }

    for ( std::size_t gid = 0; gid < globalFaceRowCount; ++ gid )
    {
        if ( ! referenced[ gid ] )
        {
            outUnreferencedSorted.push_back( static_cast< int >( gid ) );
        }
    }
}

EndNameSpace

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
    const LinkField & childNodeIds = splitTable.childFaceNodeIds;
    const LinkField & childNodeFlags = splitTable.childFaceNodeFlags;

    for ( int childIndex = 0; childIndex < childFaceIds.size(); ++ childIndex )
    {
        const int childFaceId = childFaceIds[ childIndex ];
        if ( childFaceId < 0 ||
             static_cast< size_t >( childFaceId ) >= childNodeIds.size() ||
             static_cast< size_t >( childFaceId ) >= childNodeFlags.size() )
        {
            throw std::out_of_range(
                "GetInterfaceChildFaceIds: child face ID is out of range" );
        }
        if ( childNodeIds[ childFaceId ].size() !=
             childNodeFlags[ childFaceId ].size() )
        {
            throw std::logic_error(
                "GetInterfaceChildFaceIds: child face node IDs and flags have different sizes" );
        }
    }

    return childFaceIds;
}

const IntField & GetInterfaceChildFaceIds(
    const IFaceLink & interfaceLink, int zoneId, int localInterfaceFaceId )
{
    if ( ! interfaceLink.face_search )
    {
        throw std::logic_error(
            "GetInterfaceChildFaceIds: face search is not initialized" );
    }
    return GetInterfaceChildFaceIds(
        interfaceLink.l2g,
        interfaceLink.face_search->interfaceSplit,
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

    int nBFaces = static_cast< int >( nBoundaryFaces );
    this->rebuildFaces.resize( 0 );

    int localInterfaceFaceId = 0;
    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int bcType = this->bcManager->bcRecord->bcType[ iFace ];

        int nFNode = this->faces[ iFace ].size();
        IntField resolvedFaceNodes;

        if ( BC::IsInterfaceBc( bcType ) )
        {
            const IntField & childFaceIds = GetInterfaceChildFaceIds( iFaceLink.l2g, iFaceLink.face_search->interfaceSplit, this->GetGrid().id, localInterfaceFaceId );
            int nCFace = childFaceIds.size();

            if ( nCFace > 0 )
            {
                for ( int iCFace = 0; iCFace < nCFace; ++ iCFace )
                {
                    //this->lCellNew.push_back( this->lCell[ iFace ] );

                    const int cFid = childFaceIds[ iCFace ];
                    const bool hasNodeMesh = static_cast< bool >( this->GetGrid().nodeMesh );
                    const HXSize_t nodeMeshCount = hasNodeMesh
                        ? this->GetGrid().nodeMesh->GetNumberOfNodes()
                        : 0;
                    resolvedFaceNodes = ResolveChildFaceNodeIds(
                        this->faces[ iFace ],
                        iFaceLink.face_search->interfaceSplit.childFaceNodeIds[ cFid ],
                        iFaceLink.face_search->interfaceSplit.childFaceNodeFlags[ cFid ],
                        hasNodeMesh,
                        nodeMeshCount );
                    this->rebuildFaces.push_back( resolvedFaceNodes );
                }
            }
            else
            {
                for ( int iNode = 0; iNode < nFNode; ++ iNode )
                {
                    int nodeIndex = this->faces[ iFace ][ iNode ];
                    resolvedFaceNodes.push_back( nodeIndex );
                }
                this->rebuildFaces.push_back( resolvedFaceNodes );
            }
            ++ localInterfaceFaceId;
        }
        else
        {
            for ( int iNode = 0; iNode < nFNode; ++ iNode )
            {
                int nodeIndex = this->faces[ iFace ][ iNode ];
                resolvedFaceNodes.push_back( nodeIndex );
            }
            this->rebuildFaces.push_back( resolvedFaceNodes );
        }
    }

    //Inner Face
    for ( int iFace = nBFaces; iFace < nFaces; ++ iFace )
    {
        int nFNode = this->faces[ iFace ].size();
        IntField resolvedFaceNodes;
        for ( int iNode = 0; iNode < nFNode; ++ iNode )
        {
            int nId = this->faces[ iFace ][ iNode ];
            resolvedFaceNodes.push_back( nId );
        }
        this->rebuildFaces.push_back( resolvedFaceNodes );
    }
}

void FaceTopo::SetNewFace2Cell( IFaceLink & iFaceLink )
{
    int nBFaces = this->bcManager->bcRecord->GetBoundaryFaceCount();

    int localInterfaceFaceId = 0;

    this->rebuildLeftCells.resize( 0 );
    this->rebuildRightCells.resize( 0 );

    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int bcType = this->bcManager->bcRecord->bcType[ iFace ];

        if ( BC::IsInterfaceBc( bcType ) )
        {
            const IntField & childFaceIds = GetInterfaceChildFaceIds( iFaceLink.l2g, iFaceLink.face_search->interfaceSplit, this->GetGrid().id, localInterfaceFaceId );
            int nCFace = childFaceIds.size();

            if ( nCFace > 0 )
            {
                for ( int iChildFace = 0; iChildFace < nCFace; ++ iChildFace )
                {
                    this->rebuildLeftCells.push_back( this->lCells[ iFace ] );
                }
            }
            else
            {
                this->rebuildLeftCells.push_back( this->lCells[ iFace ] );
            }
            ++ localInterfaceFaceId;
        }
        else
        {
            this->rebuildLeftCells.push_back( this->lCells[ iFace ] );
        }
    }

    int nBFaceNew = this->rebuildLeftCells.size();

    int nCells = this->GetGrid().nCells;

    for ( int iFace = 0; iFace < nBFaceNew; ++ iFace )
    {
        this->rebuildRightCells.push_back( iFace + nCells );
    }

    int nFaces = this->GetFaceCount();

    for ( int iFace = nBFaces; iFace < nFaces; ++ iFace )
    {
        this->rebuildLeftCells.push_back( this->lCells[ iFace ] );
        this->rebuildRightCells.push_back( this->rCells[ iFace ] );
    }
}


void FaceTopo::ModifyBoundaryInformation( IFaceLink & iFaceLink )
{
    //find the new number of boundary faces
    int nBFaces = this->bcManager->bcRecord->GetBoundaryFaceCount();

    int localIid = 0;

    IntField localI2B;

    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int bcType = this->bcManager->bcRecord->bcType[ iFace ];
        if ( BC::IsInterfaceBc( bcType ) )
        {
            localI2B.push_back( iFace );
            ++ localIid;
        }
    }

    int nBFaceNew = 0;

    int nIFaces = localI2B.size();

    int nIFaceNew = nIFaces;

    iFaceLink.nChild.resize( iFaceLink.l2g.size() );
    iFaceLink.nChild[ this->GetGrid().id ].assign( nIFaces, 0 );

    for ( int iFid = 0; iFid < nIFaces; ++ iFid )
    {
        const IntField & childFaceIds = GetInterfaceChildFaceIds( iFaceLink.l2g, iFaceLink.face_search->interfaceSplit, this->GetGrid().id, iFid );
        const int gFid = iFaceLink.l2g[ this->GetGrid().id ][ iFid ];
        int nCFace = childFaceIds.size();

        if ( nCFace > 0 )
        {
            iFaceLink.nChild[ this->GetGrid().id ][ iFid ] = nCFace;
            for ( int iCFace = 0; iCFace < nCFace; ++ iCFace )
            {
                int cFid = childFaceIds[ iCFace ];
                iFaceLink.l2gNew[ this->GetGrid().id ].push_back( cFid );

                iFaceLink.nChild[ this->GetGrid().id ].push_back( 0 );
            }
            ++ nIFaceNew;
        }
        else
        {
            iFaceLink.l2gNew[ this->GetGrid().id ].push_back( gFid );
        }
    }

    this->ResetNumberOfBoundaryCondition( iFaceLink );
}

void FaceTopo::ResetNumberOfBoundaryCondition( IFaceLink & iFaceLink )
{
    int nBFaces = this->bcManager->bcRecord->GetBoundaryFaceCount();

    this->bcManager->bcRecordNew->bcType.resize( 0 );
    this->bcManager->bcRecordNew->bcNameId.resize( 0 );

    int localIid = 0;
    for ( int iFace = 0; iFace < nBFaces; ++ iFace )
    {
        int bcRegion = this->bcManager->bcRecord->bcNameId[ iFace ];
        int bcType = this->bcManager->bcRecord->bcType[ iFace ];
        if ( BC::IsInterfaceBc( bcType ) )
        {
            const IntField & childFaceIds = GetInterfaceChildFaceIds( iFaceLink.l2g, iFaceLink.face_search->interfaceSplit, this->GetGrid().id, localIid );
            int nCFace = childFaceIds.size();

            if ( nCFace > 0 )
            {
                for ( int iCFace = 0; iCFace < nCFace; ++ iCFace )
                {
                    this->bcManager->bcRecordNew->bcType.push_back( bcType );
                    this->bcManager->bcRecordNew->bcNameId.push_back( bcRegion );
                }
            }
            else
            {
                this->bcManager->bcRecordNew->bcType.push_back( bcType );
                this->bcManager->bcRecordNew->bcNameId.push_back( bcRegion );
            }

            ++ localIid;
        }
        else
        {
            this->bcManager->bcRecordNew->bcType.push_back( bcType );
            this->bcManager->bcRecordNew->bcNameId.push_back( bcRegion );
        }
    }

    this->ConstructNewInterfaceMap( iFaceLink );
}

void FaceTopo::ConstructNewInterfaceMap( IFaceLink & iFaceLink )
{
    int nIFaceNew = iFaceLink.l2gNew[ this->GetGrid().id ].size();

    for ( int localIFid = 0; localIFid < nIFaceNew; ++ localIFid )
    {
        int gIFid = iFaceLink.l2gNew[ this->GetGrid().id ][ localIFid ];
        int oldSize = iFaceLink.gI2ZidNew.size();
        int newSize = ONEFLOW::MAX( gIFid + 1, oldSize );

        iFaceLink.gI2ZidNew.resize( newSize );
        iFaceLink.g2lNew.resize( newSize );

        int nFZone = iFaceLink.gI2ZidNew[ gIFid ].size();

        if ( nFZone < 2 )
        {
            iFaceLink.gI2ZidNew[ gIFid ].push_back( this->GetGrid().id );
            iFaceLink.g2lNew   [ gIFid ].push_back( localIFid );
        }
    }
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

EndNameSpace

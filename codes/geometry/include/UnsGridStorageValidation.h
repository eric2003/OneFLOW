#pragma once

#include "UnsGrid.h"
#include "FaceTopo.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include "BcRecord.h"
#include "InterFace.h"

#include <stdexcept>
#include <unordered_set>

namespace ONEFLOW
{

// Validate stored mesh counts independently from the metadata adapter.
inline void ValidateUnsGridStorageCounts( const UnsGrid & grid )
{
    const bool hasNodeStorage = static_cast< bool >( grid.nodeMesh );
    const bool hasMeshStorage = grid.HasMeshStorage();

    if ( ! hasNodeStorage && ! hasMeshStorage )
    {
        if ( grid.nNodes != 0 || grid.nFaces != 0 || grid.nCells != 0 ||
             grid.nBFaces != 0 || grid.nIFaces != 0 )
        {
            throw std::logic_error(
                "ValidateUnsGridStorageCounts: counts are nonzero after storage release" );
        }
        return;
    }

    if ( ! hasNodeStorage || ! hasMeshStorage )
    {
        throw std::logic_error(
            "ValidateUnsGridStorageCounts: grid storage is only partially initialized" );
    }

    if ( grid.nNodes < 0 || grid.nFaces < 0 || grid.nCells < 0 ||
         grid.nBFaces < 0 || grid.nIFaces < 0 )
    {
        throw std::logic_error(
            "ValidateUnsGridStorageCounts: mesh counts must be non-negative" );
    }

    if ( grid.nodeMesh->GetNumberOfNodes() != static_cast< HXSize_t >( grid.nNodes ) )
    {
        throw std::logic_error(
            "ValidateUnsGridStorageCounts: node count does not match coordinate storage" );
    }

    const FaceTopo & faceTopo = grid.GetFaceTopo();
    const CellTopo & cellTopo = grid.GetCellMesh().GetCellTopo();
    const HXSize_t faceCount = static_cast< HXSize_t >( grid.nFaces );
    const HXSize_t cellCount = static_cast< HXSize_t >( grid.nCells );
    const HXSize_t boundaryFaceCount = faceTopo.GetNBFaces();

    if ( faceTopo.GetNFaces() != faceCount ||
         faceTopo.GetFaces().size() != faceCount ||
         faceTopo.GetLeftCells().size() != faceCount ||
         faceTopo.GetRightCells().size() != faceCount )
    {
        throw std::logic_error(
            "ValidateUnsGridStorageCounts: face count does not match face topology storage" );
    }

    for ( HXSize_t faceId = 0; faceId < faceCount; ++ faceId )
    {
        for ( int nodeId : faceTopo.GetFaces()[ faceId ] )
        {
            if ( nodeId < 0 ||
                 static_cast< HXSize_t >( nodeId ) >=
                     static_cast< HXSize_t >( grid.nNodes ) )
            {
                throw std::logic_error(
                    "ValidateUnsGridStorageCounts: face references an invalid node ID" );
            }
        }
    }

    if ( cellTopo.GetNumberOfCells() != cellCount ||
         cellTopo.elements.size() != cellCount ||
         cellTopo.blank.size() != cellCount )
    {
        throw std::logic_error(
            "ValidateUnsGridStorageCounts: cell count does not match cell topology storage" );
    }

    if ( boundaryFaceCount != static_cast< HXSize_t >( grid.nBFaces ) ||
         boundaryFaceCount > faceCount )
    {
        throw std::logic_error(
            "ValidateUnsGridStorageCounts: boundary-face count does not match face topology storage" );
    }

    const BcRecord & boundaryRecord = faceTopo.GetBcRecord();
    if ( boundaryRecord.bcType.size() != boundaryFaceCount ||
         boundaryRecord.bcNameId.size() != boundaryFaceCount )
    {
        throw std::logic_error(
            "ValidateUnsGridStorageCounts: boundary record count does not match boundary-face count" );
    }

    if ( grid.nIFaces < 0 || grid.nIFaces > grid.nBFaces || ! grid.interFace ||
         grid.interFace->nIFaces != grid.nIFaces ||
         grid.interFace->zoneId.size() != static_cast< HXSize_t >( grid.nIFaces ) ||
         grid.interFace->localInterfaceId.size() != static_cast< HXSize_t >( grid.nIFaces ) ||
         grid.interFace->localCellId.size() != static_cast< HXSize_t >( grid.nIFaces ) ||
         grid.interFace->i2b.size() != static_cast< HXSize_t >( grid.nIFaces ) ||
         grid.interFace->idir.size() != static_cast< HXSize_t >( grid.nIFaces ) )
    {
        throw std::logic_error(
            "ValidateUnsGridStorageCounts: interface count does not match interface mapping storage" );
    }

    std::unordered_set< int > mappedBoundaryFaces;
    for ( int interfaceFaceId = 0; interfaceFaceId < grid.nIFaces; ++ interfaceFaceId )
    {
        const int boundaryFaceId = grid.interFace->i2b[ interfaceFaceId ];
        if ( boundaryFaceId < 0 ||
             static_cast< HXSize_t >( boundaryFaceId ) >= boundaryFaceCount ||
             ! mappedBoundaryFaces.insert( boundaryFaceId ).second )
        {
            throw std::logic_error(
                "ValidateUnsGridStorageCounts: interface-to-boundary mapping is invalid" );
        }
    }
}

} // namespace ONEFLOW

#pragma once

#include "UnsGrid.h"
#include "FaceTopo.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include "BcRecord.h"

#include <stdexcept>

namespace ONEFLOW
{

// Validate stored mesh counts independently from the metadata adapter.
inline void ValidateUnsGridStorageCounts( const UnsGrid & grid )
{
    const bool hasNodeStorage = static_cast< bool >( grid.nodeMesh );
    const bool hasFaceStorage = static_cast< bool >( grid.faceMesh );
    const bool hasCellStorage = static_cast< bool >( grid.cellMesh );

    if ( ! hasNodeStorage && ! hasFaceStorage && ! hasCellStorage )
    {
        if ( grid.nNodes != 0 || grid.nFaces != 0 || grid.nCells != 0 ||
             grid.nBFaces != 0 )
        {
            throw std::logic_error(
                "ValidateUnsGridStorageCounts: counts are nonzero after storage release" );
        }
        return;
    }

    if ( ! hasNodeStorage || ! hasFaceStorage || ! hasCellStorage )
    {
        throw std::logic_error(
            "ValidateUnsGridStorageCounts: grid storage is only partially initialized" );
    }

    if ( grid.nNodes < 0 || grid.nodeMesh->GetNumberOfNodes() !=
        static_cast< HXSize_t >( grid.nNodes ) )
    {
        throw std::logic_error(
            "ValidateUnsGridStorageCounts: node count does not match coordinate storage" );
    }

    const FaceTopo & faceTopo = grid.GetFaceTopo();
    const CellTopo & cellTopo = grid.GetCellMesh().GetCellTopo();
    const HXSize_t faceCount = static_cast< HXSize_t >( grid.nFaces );
    const HXSize_t cellCount = static_cast< HXSize_t >( grid.nCells );
    const HXSize_t boundaryFaceCount = faceTopo.GetNBFaces();

    if ( grid.nFaces < 0 || faceTopo.GetNFaces() != faceCount ||
         faceTopo.GetFaces().size() != faceCount ||
         faceTopo.GetLeftCells().size() != faceCount ||
         faceTopo.GetRightCells().size() != faceCount )
    {
        throw std::logic_error(
            "ValidateUnsGridStorageCounts: face count does not match face topology storage" );
    }

    if ( grid.nCells < 0 || cellTopo.GetNumberOfCells() != cellCount ||
         cellTopo.elements.size() != cellCount ||
         cellTopo.blank.size() != cellCount )
    {
        throw std::logic_error(
            "ValidateUnsGridStorageCounts: cell count does not match cell topology storage" );
    }

    if ( grid.nBFaces < 0 ||
         boundaryFaceCount != static_cast< HXSize_t >( grid.nBFaces ) ||
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
}

} // namespace ONEFLOW

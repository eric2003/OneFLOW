#pragma once

#include "Grid.h"
#include "GridSummary.h"
#include "NodeMesh.h"
#include "UnsGrid.h"
#include "FaceTopo.h"
#include "CellMesh.h"
#include "CellTopo.h"
#include "BcRecord.h"

#include <stdexcept>

namespace ONEFLOW
{

// Capture L0 metadata without exposing the grid storage to serializers.
inline GridSummary MakeGridSummary( const Grid & grid )
{
    GridSummary summary{
        grid.dimension,
        grid.nNodes,
        grid.nFaces,
        grid.nCells,
        grid.nBFaces
    };

    if ( grid.nodeMesh &&
         ( grid.nNodes < 0 ||
           grid.nodeMesh->GetNumberOfNodes() !=
               static_cast< HXSize_t >( grid.nNodes ) ) )
    {
        throw std::invalid_argument(
            "MakeGridSummary: node count does not match coordinate storage" );
    }

    if ( ! IsValidGridSummary( summary ) )
    {
        throw std::invalid_argument( "MakeGridSummary: grid metadata is inconsistent" );
    }

    return summary;
}

// Check storage counts separately from metadata extraction so callers can
// choose the level of validation appropriate to their workflow.
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

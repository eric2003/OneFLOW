#pragma once

#include "Grid.h"
#include "GridSummary.h"
#include "NodeMesh.h"

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

} // namespace ONEFLOW

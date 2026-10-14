#pragma once

#include "Grid.h"
#include "GridSummary.h"

namespace ONEFLOW
{

// Capture L0 metadata without exposing the grid storage to serializers.
inline GridSummary MakeGridSummary( const Grid & grid )
{
    return GridSummary{
        grid.dimension,
        grid.nNodes,
        grid.nFaces,
        grid.nCells,
        grid.nBFaces
    };
}

} // namespace ONEFLOW

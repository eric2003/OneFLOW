#include "GridInfoSnapshotAdapter.h"

#include "Grid.h"
#include "NodeMesh.h"

#include <algorithm>
#include <limits>

namespace ONEFLOW
{
GridInfoSnapshot MakeGridInfoSnapshot( const Grids & grids )
{
    GridInfoSnapshot snapshot;
    snapshot.dimension = grids.empty() ? 0 : GridAt( grids, 0 ).dimension;
    snapshot.zones.reserve( grids.size() );

    std::array< double, 3 > minimum = {
        std::numeric_limits< double >::max(),
        std::numeric_limits< double >::max(),
        std::numeric_limits< double >::max()
    };
    std::array< double, 3 > maximum = {
        std::numeric_limits< double >::lowest(),
        std::numeric_limits< double >::lowest(),
        std::numeric_limits< double >::lowest()
    };

    for ( const auto & gridPtr : grids )
    {
        const Grid & grid = *gridPtr;
        snapshot.nodeCount += grid.nNodes;
        snapshot.faceCount += grid.nFaces;
        snapshot.cellCount += grid.nCells;
        snapshot.boundaryFaceCount += grid.nBFaces;

        snapshot.zones.push_back( {
            grid.name,
            grid.dimension,
            grid.nNodes,
            grid.nFaces,
            grid.nCells,
            grid.nBFaces
        } );

        if ( ! grid.nodeMesh ) continue;
        const NodeMesh & nodes = *grid.nodeMesh;
        for ( HXSize_t iNode = 0; iNode < nodes.GetNumberOfNodes(); ++ iNode )
        {
            const NodeMesh::Point3 point = nodes.GetPoint( iNode );
            const std::array< double, 3 > coordinates = {
                static_cast< double >( point.x ),
                static_cast< double >( point.y ),
                static_cast< double >( point.z )
            };
            for ( size_t axis = 0; axis < coordinates.size(); ++ axis )
            {
                minimum[ axis ] = std::min( minimum[ axis ], coordinates[ axis ] );
                maximum[ axis ] = std::max( maximum[ axis ], coordinates[ axis ] );
            }
            snapshot.hasBoundingBox = true;
        }
    }

    if ( snapshot.hasBoundingBox )
    {
        snapshot.minimum = minimum;
        snapshot.maximum = maximum;
    }
    return snapshot;
}
} // namespace ONEFLOW

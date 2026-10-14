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
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.
\*---------------------------------------------------------------------------*/

#include "GridInfoWriter.h"

#include "Grid.h"
#include "NodeMesh.h"
#include "Prj.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <nlohmann/json.hpp>

namespace ONEFLOW
{
void WriteGridInfo( const Grids & grids, const GridConfig & config )
{
    if ( config.gridInfoFile.empty() ) return;

    std::array< Real, 3 > minimum = {
        std::numeric_limits< Real >::max(),
        std::numeric_limits< Real >::max(),
        std::numeric_limits< Real >::max()
    };
    std::array< Real, 3 > maximum = {
        std::numeric_limits< Real >::lowest(),
        std::numeric_limits< Real >::lowest(),
        std::numeric_limits< Real >::lowest()
    };
    bool hasCoordinates = false;
    int nodeCount = 0;
    int faceCount = 0;
    int cellCount = 0;
    int boundaryFaceCount = 0;
    const int dimension = grids.empty() ? 0 : GridAt( grids, 0 ).dimension;

    nlohmann::ordered_json zones = nlohmann::ordered_json::array();
    for ( const auto & gridPtr : grids )
    {
        const Grid & grid = *gridPtr;
        nodeCount += grid.nNodes;
        faceCount += grid.nFaces;
        cellCount += grid.nCells;
        boundaryFaceCount += grid.nBFaces;

        nlohmann::ordered_json zone;
        zone[ "name" ] = grid.name;
        zone[ "dimension" ] = grid.dimension;
        zone[ "nodes" ] = grid.nNodes;
        zone[ "faces" ] = grid.nFaces;
        zone[ "cells" ] = grid.nCells;
        zone[ "boundaryFaces" ] = grid.nBFaces;
        zones.push_back( std::move( zone ) );

        if ( ! grid.nodeMesh ) continue;
        const NodeMesh & nodes = *grid.nodeMesh;
        for ( HXSize_t iNode = 0; iNode < nodes.GetNumberOfNodes(); ++ iNode )
        {
            const NodeMesh::Point3 nodePoint = nodes.GetPoint( iNode );
            const std::array< Real, 3 > point = {
                nodePoint.x, nodePoint.y, nodePoint.z
            };
            for ( size_t axis = 0; axis < point.size(); ++ axis )
            {
                minimum[ axis ] = std::min( minimum[ axis ], point[ axis ] );
                maximum[ axis ] = std::max( maximum[ axis ], point[ axis ] );
            }
            hasCoordinates = true;
        }
    }

    nlohmann::ordered_json document;
    document[ "schema" ] = "oneflow.mesh-info";
    document[ "version" ] = 1;
    document[ "dimension" ] = dimension;
    document[ "counts" ] = {
        { "zones", static_cast< int >( grids.size() ) },
        { "nodes", nodeCount },
        { "faces", faceCount },
        { "cells", cellCount },
        { "boundaryFaces", boundaryFaceCount }
    };

    if ( hasCoordinates )
    {
        document[ "boundingBox" ] = {
            { "min", { minimum[ 0 ], minimum[ 1 ], minimum[ 2 ] } },
            { "max", { maximum[ 0 ], maximum[ 1 ], maximum[ 2 ] } }
        };
    }
    else
    {
        document[ "boundingBox" ] = nullptr;
    }
    document[ "zones" ] = std::move( zones );

    const std::string outputPath = Prj::GetPrjFileName( config.gridInfoFile );
    std::ofstream output( outputPath, std::ios_base::out | std::ios_base::trunc );
    output << document.dump( 2 ) << '\n';
    if ( ! output )
    {
        throw std::runtime_error(
            "WriteGridInfo: failed to write grid information file: " + outputPath );
    }
}
} // namespace ONEFLOW

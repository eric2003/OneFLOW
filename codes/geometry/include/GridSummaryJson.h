#pragma once

#include "GridSummary.h"

#include <nlohmann/json.hpp>

#include <ostream>
#include <stdexcept>

namespace ONEFLOW
{

// Writes the versioned L0 metadata document. This does not serialize mesh topology.
inline void WriteGridSummaryJson( std::ostream & output, const GridSummary & summary )
{
    const std::vector< std::string > errors = ValidateGridSummary( summary );
    if ( ! errors.empty() )
    {
        throw std::invalid_argument( "WriteGridSummaryJson: invalid grid summary" );
    }

    nlohmann::ordered_json document;
    document[ "schema" ] = "oneflow.grid-summary";
    document[ "version" ] = 1;
    document[ "dimension" ] = summary.dimension;

    nlohmann::ordered_json counts;
    counts[ "nodes" ] = summary.nodeCount;
    counts[ "faces" ] = summary.faceCount;
    counts[ "cells" ] = summary.cellCount;
    counts[ "boundaryFaces" ] = summary.boundaryFaceCount;
    document[ "counts" ] = counts;

    output << document.dump( 2 ) << '\n';
}

} // namespace ONEFLOW

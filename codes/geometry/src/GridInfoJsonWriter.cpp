#include "GridInfoJsonWriter.h"

#include <nlohmann/json.hpp>

#include <ostream>
#include <stdexcept>

namespace ONEFLOW
{
void WriteGridInfoJson( std::ostream & output, const GridInfoSnapshot & snapshot )
{
    nlohmann::ordered_json zones = nlohmann::ordered_json::array();
    for ( const GridInfoZoneSnapshot & zone : snapshot.zones )
    {
        zones.push_back( {
            { "name", zone.name },
            { "dimension", zone.dimension },
            { "nodes", zone.nodeCount },
            { "faces", zone.faceCount },
            { "cells", zone.cellCount },
            { "boundaryFaces", zone.boundaryFaceCount }
        } );
    }

    nlohmann::ordered_json document;
    document[ "schema" ] = "oneflow.mesh-info";
    document[ "version" ] = 1;
    document[ "dimension" ] = snapshot.dimension;
    document[ "counts" ] = {
        { "zones", snapshot.zones.size() },
        { "nodes", snapshot.nodeCount },
        { "faces", snapshot.faceCount },
        { "cells", snapshot.cellCount },
        { "boundaryFaces", snapshot.boundaryFaceCount }
    };

    if ( snapshot.hasBoundingBox )
    {
        document[ "boundingBox" ] = {
            { "min", snapshot.minimum },
            { "max", snapshot.maximum }
        };
    }
    else
    {
        document[ "boundingBox" ] = nullptr;
    }
    document[ "zones" ] = std::move( zones );

    output << document.dump( 2 ) << '\\n';
    if ( ! output )
    {
        throw std::runtime_error( "WriteGridInfoJson: failed to write JSON data" );
    }
}
} // namespace ONEFLOW

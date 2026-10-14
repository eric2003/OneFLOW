#pragma once

#include "GridSummary.h"

#include <ostream>

namespace ONEFLOW
{

// Writes the versioned L0 metadata document. This does not serialize mesh topology.
inline void WriteGridSummaryJson( std::ostream & output, const GridSummary & summary )
{
    output << "{\n"
           << "  \"schema\": \"oneflow.grid-summary\",\n"
           << "  \"version\": 1,\n"
           << "  \"dimension\": " << summary.dimension << ",\n"
           << "  \"counts\": {\n"
           << "    \"nodes\": " << summary.nodeCount << ",\n"
           << "    \"faces\": " << summary.faceCount << ",\n"
           << "    \"cells\": " << summary.cellCount << ",\n"
           << "    \"boundaryFaces\": " << summary.boundaryFaceCount << "\n"
           << "  }\n"
           << "}\n";
}

} // namespace ONEFLOW

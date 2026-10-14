#pragma once

#include <string>
#include <vector>

namespace ONEFLOW
{

// L0 metadata for describing a grid without owning its mesh storage.
struct GridSummary
{
    int dimension = 0;
    int nodeCount = 0;
    int faceCount = 0;
    int cellCount = 0;
    int boundaryFaceCount = 0;
};

// Returns all contract violations without requiring Grid initialization or globals.
inline std::vector< std::string > ValidateGridSummary( const GridSummary & summary )
{
    std::vector< std::string > errors;

    if ( summary.dimension < 1 || summary.dimension > 3 )
    {
        errors.emplace_back( "dimension must be between 1 and 3" );
    }
    if ( summary.nodeCount < 0 )
    {
        errors.emplace_back( "node count must not be negative" );
    }
    if ( summary.faceCount < 0 )
    {
        errors.emplace_back( "face count must not be negative" );
    }
    if ( summary.cellCount < 0 )
    {
        errors.emplace_back( "cell count must not be negative" );
    }
    if ( summary.boundaryFaceCount < 0 )
    {
        errors.emplace_back( "boundary face count must not be negative" );
    }
    if ( summary.boundaryFaceCount > summary.faceCount )
    {
        errors.emplace_back( "boundary face count must not exceed face count" );
    }
    if ( ( summary.faceCount > 0 || summary.cellCount > 0 ) &&
         summary.nodeCount == 0 )
    {
        errors.emplace_back( "faces or cells require at least one node" );
    }

    return errors;
}

inline bool IsValidGridSummary( const GridSummary & summary )
{
    return ValidateGridSummary( summary ).empty();
}

} // namespace ONEFLOW

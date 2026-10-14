#pragma once

#include "GridHandles.h"
#include "GridTypes.h"

namespace ONEFLOW
{

// Write runtime mesh metadata without mixing serialization into CalcGrid workflow.
void WriteGridInfo( const Grids & grids, const GridConfig & config );

} // namespace ONEFLOW

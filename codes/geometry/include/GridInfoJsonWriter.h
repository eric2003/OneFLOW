#pragma once

#include "GridInfoSnapshot.h"

#include <iosfwd>

namespace ONEFLOW
{

void WriteGridInfoJson( std::ostream & output, const GridInfoSnapshot & snapshot );

} // namespace ONEFLOW

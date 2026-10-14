#include "GridInfoWriter.h"

#include "GridInfoJsonWriter.h"
#include "GridInfoSnapshotAdapter.h"
#include "Prj.h"

#include <fstream>
#include <stdexcept>
#include <string>

namespace ONEFLOW
{
void WriteGridInfo( const Grids & grids, const GridConfig & config )
{
    if ( config.gridInfoFile.empty() ) return;

    const GridInfoSnapshot snapshot = MakeGridInfoSnapshot( grids );
    const std::string outputPath = Prj::GetPrjFileName( config.gridInfoFile );
    std::ofstream output( outputPath, std::ios_base::out | std::ios_base::trunc );
    if ( ! output )
    {
        throw std::runtime_error(
            "WriteGridInfo: failed to open grid information file: " + outputPath );
    }
    WriteGridInfoJson( output, snapshot );
}
} // namespace ONEFLOW

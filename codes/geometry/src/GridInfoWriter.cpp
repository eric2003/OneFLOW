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
    const std::string & summaryFile = config.output.summaryFile.empty()
        ? config.gridInfoFile
        : config.output.summaryFile;
    if ( summaryFile.empty() ) return;

    GridInfoSnapshot snapshot = MakeGridInfoSnapshot( grids );
    snapshot.caseRoot = Prj::prjBaseDir;
    snapshot.projectRoot = Prj::projectBaseDir;
    snapshot.inputRoot = Prj::GetInputBaseDir();
    snapshot.outputRoot = Prj::GetOutputBaseDir();

    if ( ! config.sourceFile.empty() )
    {
        const std::string sourceCaseDir =
            config.ResolveSourceCaseDir( Prj::prjBaseDir );
        snapshot.sourceGridPath =
            Prj::GetCaseFileName( sourceCaseDir, config.sourceFile );
    }

    const std::string outputGridFile = config.targetFile.empty()
        ? config.partitionFile
        : config.targetFile;
    if ( ! outputGridFile.empty() )
    {
        snapshot.outputGridPath = Prj::GetOutputFileName( outputGridFile );
    }

    const std::string outputPath = Prj::GetOutputFileName( summaryFile );
    std::ofstream output( outputPath, std::ios_base::out | std::ios_base::trunc );
    if ( ! output )
    {
        throw std::runtime_error(
            "WriteGridInfo: failed to open grid information file: " + outputPath );
    }
    WriteGridInfoJson( output, snapshot );
}
} // namespace ONEFLOW

#include <gtest/gtest.h>

#include "GridConversion.h"
#include "GridTypes.h"
#include "Prj.h"
#include "Simulation.h"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
class TemporaryCase
{
public:
    TemporaryCase()
    {
        const auto stamp =
            std::chrono::steady_clock::now().time_since_epoch().count();
        path_ = std::filesystem::temp_directory_path() /
            ("oneflow-gmsh-single-cell-" + std::to_string(stamp));
        std::filesystem::create_directories(path_);
    }

    ~TemporaryCase()
    {
        std::error_code error;
        std::filesystem::remove_all(path_, error);
    }

    const std::filesystem::path & Path() const { return path_; }

private:
    std::filesystem::path path_;
};

void ReplaceAssignment(
    std::string & contents,
    const std::string & key,
    const std::string & value )
{
    std::size_t lineStart = 0;
    while ( lineStart < contents.size() )
    {
        const std::size_t lineEnd = contents.find( '\n', lineStart );
        const std::size_t end =
            lineEnd == std::string::npos ? contents.size() : lineEnd;
        const std::size_t keyPosition = contents.find( key, lineStart );
        if ( keyPosition < end )
        {
            const std::size_t equals = contents.find( '=', keyPosition );
            const std::size_t semicolon =
                equals == std::string::npos ? std::string::npos
                                            : contents.find( ';', equals );
            if ( equals == std::string::npos || semicolon == std::string::npos ||
                 semicolon > end )
            {
                throw std::runtime_error(
                    "Invalid case parameter assignment: " + key );
            }
            contents.replace(
                equals + 1, semicolon - equals - 1, " " + value + " " );
            return;
        }
        if ( lineEnd == std::string::npos ) break;
        lineStart = lineEnd + 1;
    }
    throw std::runtime_error( "Missing case parameter assignment: " + key );
}

void ConfigureOneStepCase( const std::filesystem::path & caseDirectory )
{
    const std::filesystem::path sourceCase =
        std::filesystem::path( ONEFLOW_SOURCE_DIR ) /
        "test" / "cfd" / "plateuns2dslau2";
    std::filesystem::copy(
        sourceCase / "script",
        caseDirectory / "script",
        std::filesystem::copy_options::recursive );

    std::filesystem::create_directories( caseDirectory / "grid" );
    std::filesystem::create_directories( caseDirectory / "results" );
    std::filesystem::create_directories( caseDirectory / "restart" );
    std::filesystem::create_directories( caseDirectory / "visual" );

    std::filesystem::copy_file(
        sourceCase / "grid" / "ns_bc.txt",
        caseDirectory / "grid" / "ns_bc.txt" );

    const std::filesystem::path cfdFile =
        caseDirectory / "script" / "cfd.txt";
    std::ifstream input( cfdFile );
    if ( ! input ) throw std::runtime_error( "Cannot read copied cfd.txt" );
    std::string cfd{
        std::istreambuf_iterator<char>( input ),
        std::istreambuf_iterator<char>() };
    input.close();
    ReplaceAssignment( cfd, "maxSteps", "1" );
    ReplaceAssignment( cfd, "nFieldSave", "1" );
    ReplaceAssignment( cfd, "nForceSave", "1" );
    ReplaceAssignment( cfd, "nResSave", "1" );
    ReplaceAssignment( cfd, "nVisualSave", "1" );
    ReplaceAssignment( cfd, "gridFileName", "\"grid/single_quad.oneflow.grid\"" );
    ReplaceAssignment( cfd, "idump", "0" );
    ReplaceAssignment( cfd, "showfield", "0" );

    std::ofstream output( cfdFile, std::ios::trunc );
    if ( ! output ) throw std::runtime_error( "Cannot update copied cfd.txt" );
    output << cfd;
    if ( ! output ) throw std::runtime_error( "Failed to write copied cfd.txt" );
}
} // namespace

TEST( GmshSingleCellRuntimeIntegration,
      ConvertsReadsComputesAndWritesFlowField )
{
#ifndef ENABLE_CGNS
    GTEST_SKIP() << "CGNS is not enabled in this build.";
#else
    TemporaryCase testCase;
    const std::filesystem::path caseDirectory = testCase.Path();
    ConfigureOneStepCase( caseDirectory );

    const std::filesystem::path sourceMesh =
        std::filesystem::path( ONEFLOW_SOURCE_DIR ) /
        "grid" / "gmsh_to_oneflow" / "grid" / "single_quad.msh";
    ASSERT_TRUE( std::filesystem::is_regular_file( sourceMesh ) );
    const std::filesystem::path inputMesh =
        caseDirectory / "grid" / "single_quad.msh";
    std::filesystem::copy_file( sourceMesh, inputMesh );

    // The converter and the production CLI share the same case-root contract.
    ONEFLOW::Prj::Init();
    ONEFLOW::Prj::SetPrjBaseDir( caseDirectory.string() );

    ONEFLOW::GridConfig gridConfig;
    gridConfig.sourceType = ONEFLOW::GridFileType::Gmsh;
    gridConfig.targetType = ONEFLOW::GridFileType::OneFLOW;
    gridConfig.topology = ONEFLOW::GridTopology::Unstructured;
    gridConfig.sourceFile = "grid/single_quad.msh";
    gridConfig.targetFile = "grid/single_quad.oneflow.grid";

    ASSERT_NO_THROW(
        ONEFLOW::ConvertGrid( gridConfig, caseDirectory.string() ) );
    ASSERT_TRUE( std::filesystem::is_regular_file(
        caseDirectory / "grid" / "single_quad.oneflow.grid" ) );
    ASSERT_TRUE( std::filesystem::is_regular_file(
        caseDirectory / "grid" / "bcRegionMap.txt" ) );

    // Exercise the same production command-line route as:
    // OneFLOW.exe r <case-directory>
    std::vector<std::string> arguments{
        "OneFLOW", "r", caseDirectory.string()
    };
    std::vector<char *> argv;
    argv.reserve( arguments.size() );
    for ( std::string & argument : arguments )
    {
        argv.push_back( argument.data() );
    }

    ONEFLOW::Simulation simulation(
        static_cast<int>( argv.size() ), argv.data() );
    ASSERT_EQ( simulation.Run(), 0 )
        << "OneFLOW failed while reading or advancing the one-cell Gmsh case.";

    const std::vector<std::filesystem::path> expectedOutputs{
        caseDirectory / "results" / "res.dat",
        caseDirectory / "restart" / "flow.dat",
        caseDirectory / "visual" / "visual.dat"
    };
    for ( const auto & outputPath : expectedOutputs )
    {
        ASSERT_TRUE( std::filesystem::is_regular_file( outputPath ) )
            << "Missing runtime output: " << outputPath.string();
        EXPECT_GT( std::filesystem::file_size( outputPath ), 0u )
            << "Empty runtime output: " << outputPath.string();
    }

    std::ifstream visualInput( expectedOutputs.back() );
    const std::string visual{
        std::istreambuf_iterator<char>( visualInput ),
        std::istreambuf_iterator<char>() };
    EXPECT_NE( visual.find( "Elements = 1" ), std::string::npos );
#endif
}

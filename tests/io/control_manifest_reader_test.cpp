#include "gtest/gtest.h"
#include "ParaFile.h"
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace ONEFLOW;

TEST( ControlManifestReader, PreservesListedScriptOrderAndCaseRoot )
{
    const std::filesystem::path caseDir =
        std::filesystem::temp_directory_path()
        / "oneflow_control_manifest_reader_test";
    const std::filesystem::path scriptDir = caseDir / "script";

    std::filesystem::remove_all( caseDir );
    std::filesystem::create_directories( scriptDir );

    {
        std::ofstream manifest( scriptDir / "control.txt" );
        ASSERT_TRUE( manifest.good() );
        manifest
            << "\"simu.txt\";\n"
            << "\n"
            << "\"dimension.txt\";\n"
            << "\"grid.txt\";\n";
    }

    std::vector< std::string > scriptFiles = { "existing-entry.txt" };
    ReadScriptFileNameList( caseDir.string(), scriptFiles );

    ASSERT_EQ( scriptFiles.size(), static_cast< std::size_t >( 4 ) );
    EXPECT_EQ( scriptFiles[ 0 ], "existing-entry.txt" );
    EXPECT_EQ(
        std::filesystem::path( scriptFiles[ 1 ] ),
        ( scriptDir / "simu.txt" ).lexically_normal() );
    EXPECT_EQ(
        std::filesystem::path( scriptFiles[ 2 ] ),
        ( scriptDir / "dimension.txt" ).lexically_normal() );
    EXPECT_EQ(
        std::filesystem::path( scriptFiles[ 3 ] ),
        ( scriptDir / "grid.txt" ).lexically_normal() );

    std::filesystem::remove_all( caseDir );
}

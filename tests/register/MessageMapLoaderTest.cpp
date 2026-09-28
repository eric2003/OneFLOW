

// MessageMapLoaderTest.cpp
#include <gtest/gtest.h>
#include "Prj.h"
#include "MessageMapLoader.h"
#include <filesystem>
#include <fstream>

using ONEFLOW::Prj;

// GetMsgFileNameList reads the manifest file
// (Prj::GetSystemFileName("action/actionFileList.txt")) and expands each
// listed entry into a full path under the same "action/" directory.
// Prj::system_root is a plain writable static std::string, so it can be
// redirected to a temp directory for the duration of this test.
class MessageMapLoaderTest : public ::testing::Test
{
protected:
    std::string savedRoot;
    std::filesystem::path tempRoot;

    void SetUp() override
    {
        savedRoot = Prj::system_root;

        tempRoot = std::filesystem::temp_directory_path()
            / "OneFLOW_MessageMapLoaderTest" / "system";

        std::filesystem::remove_all( tempRoot );
        std::filesystem::create_directories( tempRoot / "action" );

        // Prj::system_root keeps a trailing slash by convention (see Prj.cpp).
        Prj::system_root = tempRoot.string() + "/";
    }

    void TearDown() override
    {
        Prj::system_root = savedRoot;
        std::filesystem::remove_all( tempRoot );
    }
};

TEST_F( MessageMapLoaderTest, ExpandsManifestEntriesToFullPaths )
{
    const std::filesystem::path manifestFile =
        tempRoot / "action" / "actionFileList.txt";

    // Mix in a comment line and a blank line to exercise the
    // ReadNextMeaningfulLine() skip-blank/skip-comment behaviour.
    std::ofstream manifest( manifestFile );
    manifest << "# a comment line, must be skipped\n";
    manifest << "\n";
    manifest << "gridAction.txt\n";
    manifest << "solverAction.txt\n";
    manifest.close();

    ONEFLOW::StringField fileNameList;
    ONEFLOW::GetMsgFileNameList( fileNameList );

    ASSERT_EQ( fileNameList.size(), 2 );

    const std::filesystem::path expectedFirst =
        tempRoot / "action" / "gridAction.txt";
    const std::filesystem::path expectedSecond =
        tempRoot / "action" / "solverAction.txt";

    EXPECT_EQ(
        std::filesystem::path( fileNameList[ 0 ] ).lexically_normal(),
        expectedFirst.lexically_normal() );
    EXPECT_EQ(
        std::filesystem::path( fileNameList[ 1 ] ).lexically_normal(),
        expectedSecond.lexically_normal() );
}
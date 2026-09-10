// MessageMapLoaderTest.cpp
#include <gtest/gtest.h>
#include "Prj.h"
#include "MessageMapLoader.h"

// GetMsgFileNameList reads Prj::system_root + "action/actionFileList.txt".
// We don't have Prj's source, but if Prj::system_root is a writable
// static (not a compile-time constant), we can point it at a temp
// directory for the duration of this test and restore it afterward.
// If Prj::system_root turns out to be read-only/fixed at build time,
// this test should be deleted and this function left uncovered until
// the golden-trace integration stage, which will set up a real minimal
// project directory anyway.

class MessageMapLoaderTest : public ::testing::Test
{
protected:
    std::string savedRoot;

    void SetUp() override
    {
        savedRoot = ONEFLOW::Prj::system_root;
    }

    void TearDown() override
    {
        ONEFLOW::Prj::system_root = savedRoot;
    }
};

// NOTE: this test is left as a sketch, not verified to compile - it
// depends on Prj::system_root being reassignable and on being able to
// create the expected directory layout under a temp path. Please confirm
// Prj::system_root's actual declaration before enabling this test.
TEST_F( MessageMapLoaderTest, DISABLED_ExpandsManifestEntriesToFullPaths )
{
    // ... requires a real temp "action/actionFileList.txt" on disk;
    // left DISABLED_ until Prj::system_root's mutability is confirmed.
}
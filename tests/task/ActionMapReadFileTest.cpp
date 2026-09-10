// ActionMapReadFileTest.cpp
#include <gtest/gtest.h>
#include <fstream>
#include <cstdio>
#include "ActionMap.h"

// ActionMapImp::ReadFile depends on TextFileParser, which in turn delegates
// file opening to Prj::OpenFile (source not available to us). This makes
// these tests small integration tests rather than pure unit tests: if they
// fail on file-open, the cause may be a working-directory / Prj::OpenFile
// path resolution issue rather than a bug in ActionMapImp itself.

class ActionMapReadFileTest : public ::testing::Test
{
protected:
    std::string tempFilePath = "test_action_map_temp.txt";

    void WriteTempFile( const std::string & content )
    {
        std::ofstream ofs( tempFilePath );
        ofs << content;
        ofs.close();
    }

    void TearDown() override
    {
        std::remove( tempFilePath.c_str() );
    }
};

TEST_F( ActionMapReadFileTest, RegistersOneActionPerNonEmptyLine )
{
    WriteTempFile(
        "ComputeFlux\n"
        "UpdateBoundary\n"
        "WriteOutput\n"
    );

    ONEFLOW::ActionMapImp imp;
    imp.ReadFile( tempFilePath );

    EXPECT_EQ( imp.GetActionId( "ComputeFlux" ), 0 );
    EXPECT_EQ( imp.GetActionId( "UpdateBoundary" ), 1 );
    EXPECT_EQ( imp.GetActionId( "WriteOutput" ), 2 );
}

TEST_F( ActionMapReadFileTest, IgnoresBlankLines )
{
    WriteTempFile(
        "ComputeFlux\n"
        "\n"
        "\n"
        "UpdateBoundary\n"
    );

    ONEFLOW::ActionMapImp imp;
    imp.ReadFile( tempFilePath );

    EXPECT_EQ( imp.GetActionId( "ComputeFlux" ), 0 );
    EXPECT_EQ( imp.GetActionId( "UpdateBoundary" ), 1 );
}

TEST_F( ActionMapReadFileTest, OnlyTakesFirstWordOnALineWithExtraContent )
{
    // A line like "solve = 1 # comment" should only register "solve",
    // since space is a separator char and ReadNextWord() only consumes
    // the leading token before the first separator.
    WriteTempFile(
        "solve = 1\n"
        "output = 2\n"
    );

    ONEFLOW::ActionMapImp imp;
    imp.ReadFile( tempFilePath );

    EXPECT_EQ( imp.GetActionId( "solve" ), 0 );
    EXPECT_EQ( imp.GetActionId( "output" ), 1 );
    EXPECT_EQ( imp.GetActionId( "=" ), -1 );
    EXPECT_EQ( imp.GetActionId( "1" ), -1 );
}

TEST_F( ActionMapReadFileTest, DuplicateNameKeepsFirstRegisteredId )
{
    WriteTempFile(
        "solve\n"
        "output\n"
        "solve\n"
    );

    ONEFLOW::ActionMapImp imp;
    imp.ReadFile( tempFilePath );

    EXPECT_EQ( imp.GetActionId( "solve" ), 0 );
    EXPECT_EQ( imp.GetActionId( "output" ), 1 );
}

//// --- Characterization test for an existing quirk, NOT a correctness claim ---
//// ReadFile uses ReadNextNonEmptyLine(), which only skips BLANK lines.
//// It does NOT skip comment lines (lines starting with '#' or '//'), unlike
//// TextFileParser::ReadNextMeaningfulLine(). Since '#' is itself a separator
//// character, a comment-only line causes ReadNextWord() to return an empty
//// string, which then gets registered as a valid (if odd) action name "".
//// This test documents the CURRENT behavior so any future change to it is a
//// deliberate, visible decision rather than an accidental regression.
//TEST_F( ActionMapReadFileTest, CommentOnlyLineIsNotSkippedAndRegistersEmptyName )
//{
//    WriteTempFile(
//        "# this is a comment line\n"
//        "solve\n"
//    );
//
//    ONEFLOW::ActionMapImp imp;
//    imp.ReadFile( tempFilePath );
//
//    // Documenting current behavior: empty string got registered as id 0.
//    EXPECT_EQ( imp.GetActionId( "" ), 0 );
//    EXPECT_EQ( imp.GetActionId( "solve" ), 1 );
//}

//// --- Isolate finding #1: trailing empty-string registration at EOF ---
//// ReadFile's loop checks ReachTheEndOfFile() (i.e. std::fstream::eof())
//// BEFORE reading, but eofbit is only set AFTER a read attempt goes past
//// the last line. This means the loop runs one extra time past the last
//// real line, and that extra iteration registers an empty string "".
//// This happens on every file read, independent of comments.
//TEST_F( ActionMapReadFileTest, TrailingEmptyStringIsRegisteredAfterLastLine )
//{
//    WriteTempFile(
//        "solve\n"
//        "output\n"
//    );
//
//    ONEFLOW::ActionMapImp imp;
//    imp.ReadFile( tempFilePath );
//
//    EXPECT_EQ( imp.GetActionId( "solve" ), 0 );
//    EXPECT_EQ( imp.GetActionId( "output" ), 1 );
//    // Documenting the quirk: an extra "" is registered as id 2.
//    EXPECT_EQ( imp.GetActionId( "" ), 2 );
//}

//// --- Isolate finding #2: comment lines are NOT skipped and are parsed as data ---
//// ReadFile uses ReadNextNonEmptyLine() (blank-line skip only), not
//// ReadNextMeaningfulLine() (which also skips '#'/'//' comment lines).
//// Since '#' is a separator character, the first token AFTER the '#' is
//// what gets extracted and registered ¡ª not skipped, not empty.
//// NOTE: the exact token depends on Word::FindNextWord's leading-separator
//// handling; "this" is our best inference from the observed id ordering.
//// If this assertion fails, the failure message will reveal the true value.
//TEST_F( ActionMapReadFileTest, CommentLineIsNotSkippedAndIsParsedAsAToken )
//{
//    WriteTempFile(
//        "# this is a comment line\n"
//        "solve\n"
//    );
//
//    ONEFLOW::ActionMapImp imp;
//    imp.ReadFile( tempFilePath );
//
//    EXPECT_EQ( imp.GetActionId( "this" ), 0 );   // best inference; may need correction
//    EXPECT_EQ( imp.GetActionId( "solve" ), 1 );
//    EXPECT_EQ( imp.GetActionId( "" ), 2 );        // trailing EOF artifact, per finding #1
//}

TEST_F( ActionMapReadFileTest, CommentLineIsSkippedAfterFix )
{
    WriteTempFile(
        "# this is a comment line\n"
        "solve\n"
    );

    ONEFLOW::ActionMapImp imp;
    imp.ReadFile( tempFilePath );

    // The comment line should now be skipped entirely.
    EXPECT_EQ( imp.GetActionId( "solve" ), 0 );
    EXPECT_EQ( imp.GetActionId( "this" ), -1 );   // "this" must NOT be registered anymore
    EXPECT_EQ( imp.GetActionId( "" ), -1 );        // no trailing empty-string artifact either
}

TEST_F( ActionMapReadFileTest, NoTrailingEmptyStringAfterFix )
{
    WriteTempFile(
        "solve\n"
        "output\n"
    );

    ONEFLOW::ActionMapImp imp;
    imp.ReadFile( tempFilePath );

    EXPECT_EQ( imp.GetActionId( "solve" ), 0 );
    EXPECT_EQ( imp.GetActionId( "output" ), 1 );
    EXPECT_EQ( imp.GetActionId( "" ), -1 );   // the EOF-trailing artifact should be gone now
}
#include <gtest/gtest.h>
#include <fstream>
#include <cstdio>
#include "Message.h"

// MessageMapImp mirrors ActionMapImp's structure and past bug history,
// so these tests mirror ActionMapImpTest/ActionMapReadFileTest directly.

TEST( MessageMapImpTest, RegisterAssignsSequentialIds )
{
    ONEFLOW::MessageMapImp imp;
    imp.Register( "SolveNS" );
    imp.Register( "UpdateInterface" );

    EXPECT_EQ( imp.GetMsgId( "SolveNS" ), 0 );
    EXPECT_EQ( imp.GetMsgId( "UpdateInterface" ), 1 );
    EXPECT_EQ( imp.GetMsgName( 0 ), "SolveNS" );
}

TEST( MessageMapImpTest, UnknownNameReturnsNegativeOne )
{
    ONEFLOW::MessageMapImp imp;
    EXPECT_EQ( imp.GetMsgId( "not_registered" ), -1 );
}

TEST( MessageMapImpTest, RegisterIgnoresEmptyName )
{
    ONEFLOW::MessageMapImp imp;
    imp.Register( "" );
    EXPECT_EQ( imp.GetMsgId( "" ), -1 );
}

TEST( MessageMapImpTest, GetMsgNameWithOutOfRangeIdReturnsEmptyString )
{
    ONEFLOW::MessageMapImp imp;
    imp.Register( "SolveNS" );

    EXPECT_EQ( imp.GetMsgName( -1 ), "" );
    EXPECT_EQ( imp.GetMsgName( 999 ), "" );
}

TEST( MessageMapImpTest, ClearResetsAllState )
{
    ONEFLOW::MessageMapImp imp;
    imp.Register( "SolveNS" );
    imp.Clear();

    EXPECT_EQ( imp.GetMsgId( "SolveNS" ), -1 );

    imp.Register( "restart" );
    EXPECT_EQ( imp.GetMsgId( "restart" ), 0 );
}

class MessageMapReadFileTest : public ::testing::Test
{
protected:
    std::string tempFilePath = "test_message_map_temp.txt";

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

TEST_F( MessageMapReadFileTest, RegistersOneMessagePerNonEmptyLine )
{
    WriteTempFile(
        "SolveNS\n"
        "UpdateInterface\n"
    );

    ONEFLOW::MessageMapImp imp;
    imp.ReadFile( tempFilePath );

    EXPECT_EQ( imp.GetMsgId( "SolveNS" ), 0 );
    EXPECT_EQ( imp.GetMsgId( "UpdateInterface" ), 1 );
}

TEST_F( MessageMapReadFileTest, CommentLineIsSkipped )
{
    WriteTempFile(
        "# this is a comment line\n"
        "SolveNS\n"
    );

    ONEFLOW::MessageMapImp imp;
    imp.ReadFile( tempFilePath );

    EXPECT_EQ( imp.GetMsgId( "SolveNS" ), 0 );
    EXPECT_EQ( imp.GetMsgId( "this" ), -1 );
}

TEST_F( MessageMapReadFileTest, NoTrailingEmptyStringAfterEOF )
{
    WriteTempFile(
        "SolveNS\n"
        "UpdateInterface\n"
    );

    ONEFLOW::MessageMapImp imp;
    imp.ReadFile( tempFilePath );

    EXPECT_EQ( imp.GetMsgId( "" ), -1 );
}

// Facade-level test: confirms Init()/Free() reset semantics work the
// same way as ActionMap's, through the static MessageMap entry points.
class MessageMapFacadeTest : public ::testing::Test
{
protected:
    void SetUp() override { ONEFLOW::MessageMap::Init(); }
    void TearDown() override { ONEFLOW::MessageMap::Free(); }
};

TEST_F( MessageMapFacadeTest, FreeThenInitAgainStartsFromCleanState )
{
    ONEFLOW::MessageMap::Register( "SolveNS" );
    ONEFLOW::MessageMap::Free();
    ONEFLOW::MessageMap::Init();

    EXPECT_EQ( ONEFLOW::MessageMap::GetMsgId( "SolveNS" ), -1 );

    ONEFLOW::MessageMap::Register( "restart" );
    EXPECT_EQ( ONEFLOW::MessageMap::GetMsgId( "restart" ), 0 );
}

// ---------------------------------------------------------------------------
// Contract: name used by FieldSimuInitFlowField -> MultiSolverMultiGridTask.
// CmxTask resolves the string via MessageMap::GetMsgId before GenerateCmdList.
// No CFD; pure name<->id mapping (same as production lookup shape).
// ---------------------------------------------------------------------------
TEST( MessageMapImpTest, InitFlowFieldNameRoundTripsLikeCmxTaskLookup )
{
    // Must stay identical to kInitFlowFieldTaskName in FieldSimu.cpp
    const std::string kInitFlowFieldTaskName = "INIT_FLOWFIELD";

    ONEFLOW::MessageMapImp imp;
    imp.Register( kInitFlowFieldTaskName );

    const int id = imp.GetMsgId( kInitFlowFieldTaskName );
    EXPECT_GE( id, 0 );
    EXPECT_EQ( imp.GetMsgName( id ), kInitFlowFieldTaskName );

    // Same shape as CmxTask: unknown name must not look like a valid op
    EXPECT_EQ( imp.GetMsgId( "INIT_FLOWFIELD_TYPO" ), -1 );
}

TEST( MessageMapImpTest, PostProcessNameRoundTripsLikeCmxTaskLookup )
{
    // Must stay identical to kPostProcessTaskName in Multigrid.cpp
    const std::string kPostProcessTaskName = "POST_PROCESS";

    ONEFLOW::MessageMapImp imp;
    imp.Register( kPostProcessTaskName );

    const int id = imp.GetMsgId( kPostProcessTaskName );
    EXPECT_GE( id, 0 );
    EXPECT_EQ( imp.GetMsgName( id ), kPostProcessTaskName );
}


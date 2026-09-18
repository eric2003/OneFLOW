#include <gtest/gtest.h>
#include <fstream>
#include <cstdio>
#include <vector>
#include <string>
#include "Message.h"
#include "MessageMapLoader.h"
#include "CmxTaskNames.h"

// MessageMapImp mirrors ActionMapImp's structure and past bug history,

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
// Contract: CmxTaskNames constants <-> MessageMap name/id (same lookup shape
// as MultiSolverMultiGridTask / AddCmdToList production path).
// No CFD; pure mapping. Constants come from CmxTaskNames.h (single source).
// ---------------------------------------------------------------------------
TEST( MessageMapImpTest, InitFlowFieldNameRoundTripsLikeCmxTaskLookup )
{
    using ONEFLOW::kInitFlowFieldTaskName;

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
    using ONEFLOW::kPostProcessTaskName;

    ONEFLOW::MessageMapImp imp;
    imp.Register( kPostProcessTaskName );

    const int id = imp.GetMsgId( kPostProcessTaskName );
    EXPECT_GE( id, 0 );
    EXPECT_EQ( imp.GetMsgName( id ), kPostProcessTaskName );
}

TEST( MessageMapImpTest, CmxTaskNames_AddCmdPathRoundTrips )
{
    // Names used by RestartTaskReg / SolverImp / Ns / INs / Turb AddCmdToList.
    const std::vector<const char*> names = {
        ONEFLOW::kInitFirstTaskName,
        ONEFLOW::kInitRestartTaskName,
        ONEFLOW::kReadRestartTaskName,
        ONEFLOW::kInitInsRestartTaskName,
        ONEFLOW::kReadInsRestartTaskName,
        ONEFLOW::kInitFinalTaskName,
        ONEFLOW::kUploadInterfaceDataTaskName,
        ONEFLOW::kUpdateInterfaceDataTaskName,
        ONEFLOW::kDownloadInterfaceDataTaskName,
        ONEFLOW::kDumpResidualTaskName,
        ONEFLOW::kDumpAerodynamicTaskName,
        ONEFLOW::kDumpPressureCoeffTaskName,
        ONEFLOW::kDumpHeatfluxCoeffTaskName,
        ONEFLOW::kDumpRestartTaskName,
        ONEFLOW::kDumpLaminarPlateTaskName,
        ONEFLOW::kDumpTurbPlateTaskName,
        ONEFLOW::kVisualizationTaskName,
        ONEFLOW::kUpdateUnsteadyFlowTaskName,
    };

    ONEFLOW::MessageMapImp imp;
    for ( const char* name : names )
    {
        imp.Register( name );
    }

    for ( const char* name : names )
    {
        const int id = imp.GetMsgId( name );
        EXPECT_GE( id, 0 ) << name;
        EXPECT_EQ( imp.GetMsgName( id ), name ) << name;
    }
}

TEST( MessageMapImpTest, CmxTaskNames_TimeIntegralPathRoundTrips )
{
    const std::vector<const char*> names = {
        ONEFLOW::kCalcTimeStepTaskName,
        ONEFLOW::kCalcLhsTaskName,
        ONEFLOW::kUpdateFlowFieldTaskName,
        ONEFLOW::kCalcBoundaryTaskName,
        ONEFLOW::kZeroDqFieldTaskName,
        ONEFLOW::kInitLusgsTaskName,
        ONEFLOW::kLusgsLowerSweepTaskName,
        ONEFLOW::kExchangeInterfaceDqTaskName,
        ONEFLOW::kLusgsUpperSweepTaskName,
        ONEFLOW::kUpdateFlowFieldLusgsTaskName,
    };

    ONEFLOW::MessageMapImp imp;
    for ( const char* name : names )
    {
        imp.Register( name );
    }

    for ( const char* name : names )
    {
        const int id = imp.GetMsgId( name );
        EXPECT_GE( id, 0 ) << name;
        EXPECT_EQ( imp.GetMsgName( id ), name ) << name;
    }
}

TEST( MessageMapFacade, ContainsReflectsRegister )
{
    ONEFLOW::MessageMap::Init();
    EXPECT_FALSE( ONEFLOW::MessageMap::Contains( "UNITTEST_NO_SUCH_MSG" ) );
    ONEFLOW::MessageMap::Register( "UNITTEST_NO_SUCH_MSG" );
    EXPECT_TRUE( ONEFLOW::MessageMap::Contains( "UNITTEST_NO_SUCH_MSG" ) );
    ONEFLOW::MessageMap::Free();
}

TEST( CmxTaskNameValidation, CollectMissingWhenMapEmpty )
{
    ONEFLOW::MessageMap::Init();
    const ONEFLOW::StringField missing = ONEFLOW::CollectMissingCmxTaskNames();
    EXPECT_FALSE( missing.empty() );
    ONEFLOW::MessageMap::Free();
}

TEST( CmxTaskNameValidation, CollectMissingEmptyAfterRegisteringAll )
{
    ONEFLOW::MessageMap::Init();
    // Register whatever is missing until the set is complete.
    for ( int pass = 0; pass < 2; ++ pass )
    {
        const ONEFLOW::StringField missing = ONEFLOW::CollectMissingCmxTaskNames();
        for ( std::size_t i = 0; i < missing.size(); ++ i )
        {
            ONEFLOW::MessageMap::Register( missing[ i ] );
        }
    }
    EXPECT_TRUE( ONEFLOW::CollectMissingCmxTaskNames().empty() );
    ONEFLOW::MessageMap::Free();
}


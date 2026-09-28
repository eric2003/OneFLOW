#include <gtest/gtest.h>
#include "SimuTask.h"
#include "SimuContext.h"
#include "SimuTaskRequire.h"
#include <algorithm>
#include <string>
#include <vector>
#include <stdexcept>
#include <memory>

using namespace ONEFLOW;

namespace {

class MockTask : public ISimuTask
{
public:
    explicit MockTask( bool needsMap = false )
        : needsMap_( needsMap )
    {
    }

    void Execute( SimuContext& /*ctx*/ ) override { executed_ = true; }
    bool NeedsSystemMap() const override { return needsMap_; }
    bool Executed() const { return executed_; }

private:
    bool needsMap_ = false;
    bool executed_ = false;
};

} // namespace

TEST( TaskRegistryTest, RegisterAndCreate )
{
    auto& reg = TaskRegistry::Instance();

    const std::string name = "UnitTest_MockSolve";
    reg.Register( name, []() {
        return std::make_unique<MockTask>( true );
    } );

    EXPECT_TRUE( reg.Contains( name ) );

    auto task = reg.Create( name );
    ASSERT_NE( task, nullptr );
    EXPECT_TRUE( task->NeedsSystemMap() );

    SimuContext ctx( std::vector<std::string>{} );
    task->Execute( ctx );
}

TEST( TaskRegistryTest, CreateUnknownReturnsNull )
{
    auto task = TaskRegistry::Instance().Create( "DefinitelyNotRegistered_XYZ" );
    EXPECT_EQ( task, nullptr );
    EXPECT_FALSE( TaskRegistry::Instance().Contains( "DefinitelyNotRegistered_XYZ" ) );
}

TEST( TaskRegistryTest, GetAllRegisteredNamesContainsMock )
{
    auto& reg = TaskRegistry::Instance();
    const std::string name = "UnitTest_ListMe";
    reg.Register( name, []() {
        return std::make_unique<MockTask>();
    } );

    auto names = reg.GetAllRegisteredNames();
    EXPECT_NE( std::find( names.begin(), names.end(), name ), names.end() );
}

TEST( TaskRegistryTest, TaskEnumToStringMapping )
{
    EXPECT_EQ( TaskEnumToString( TaskEnum::SOLVE_FIELD ),      "Solve" );
    EXPECT_EQ( TaskEnumToString( TaskEnum::CREATE_GRID ),      "Grid" );
    EXPECT_EQ( TaskEnumToString( TaskEnum::CREATE_WALL_DIST ), "WallDist" );
    EXPECT_EQ( TaskEnumToString( TaskEnum::PARTITION_GRID ),   "Partition" );
    EXPECT_EQ( TaskEnumToString( TaskEnum::FUNCTION_TEST ),    "FunctionTest" );
    EXPECT_EQ( TaskEnumToString( TaskEnum::SOLVE_THEORY ),     "Theory" );
    EXPECT_EQ( TaskEnumToString( TaskEnum::TOY_MODEL ),        "ToyModel" );
    EXPECT_EQ( TaskEnumToString( TaskEnum::POST_TASK ),        "PostTask" );
}

TEST( TaskRegistryTest, CreateByTaskEnumUsesSameKeys )
{
    auto& reg = TaskRegistry::Instance();
    reg.Register( "Solve", []() {
        return std::make_unique<MockTask>( true );
    } );

    auto task = reg.Create( TaskEnum::SOLVE_FIELD );
    ASSERT_NE( task, nullptr );
    EXPECT_TRUE( task->NeedsSystemMap() );
}

// ---------------------------------------------------------------------------
// Contract for the production Solve path (mirrored by SolveFieldTask).
// Does NOT link SimuTaskReg / FieldSimu ¡ª only the registry + context contract.
// ---------------------------------------------------------------------------

namespace {

    // Same shape as production SolveFieldTask, without calling FieldSimu().
    class SolveFieldContractTask : public ISimuTask
    {
    public:
        bool NeedsSystemMap() const override { return true; }

        void Execute( SimuContext& ctx ) override
        {
            RequireSolveFieldContext( ctx );
            executed_ = true;  // ²»µ÷ÓÃ FieldSimu
        }

        bool executed_ = false;
    };

} // namespace

TEST( TaskRegistryTest, SolveContract_NeedsSystemMapAndAcceptsReadyCtx )
{
    auto& reg = TaskRegistry::Instance();
    const std::string name = "UnitTest_SolveContract";
    reg.Register( name, []() {
        return std::make_unique<SolveFieldContractTask>();
        } );

    auto task = reg.Create( name );
    ASSERT_NE( task, nullptr );
    EXPECT_TRUE( task->NeedsSystemMap() );

    SimuContext ctx( std::vector<std::string>{ "OneFLOW", "d", "test/plateuns2dslau2" } );
    ctx.MarkEnvironmentReady( true );
    ctx.SetTaskByName( "Solve" );

    EXPECT_NO_THROW( task->Execute( ctx ) );
}

TEST( TaskRegistryTest, SolveContract_RejectsEnvironmentNotReady )
{
    auto& reg = TaskRegistry::Instance();
    const std::string name = "UnitTest_SolveNotReady";
    reg.Register( name, []() {
        return std::make_unique<SolveFieldContractTask>();
        } );

    auto task = reg.Create( name );
    ASSERT_NE( task, nullptr );

    SimuContext ctx( std::vector<std::string>{} );
    // envReady_ default false
    ctx.SetTaskByName( "Solve" );

    EXPECT_THROW( task->Execute( ctx ), std::runtime_error );
}

TEST( TaskRegistryTest, SolveContract_RejectsWrongTaskName )
{
    auto& reg = TaskRegistry::Instance();
    const std::string name = "UnitTest_SolveWrongName";
    reg.Register( name, []() {
        return std::make_unique<SolveFieldContractTask>();
        } );

    auto task = reg.Create( name );
    ASSERT_NE( task, nullptr );

    SimuContext ctx( std::vector<std::string>{} );
    ctx.MarkEnvironmentReady( true );
    ctx.SetTaskByName( "Grid" );  // not Solve

    EXPECT_THROW( task->Execute( ctx ), std::runtime_error );
}

TEST( TaskRegistryTest, SolveEnumMapsToSolveString )
{
    // Documents the control-file key used by plateuns2dslau2 (simutask = "Solve").
    EXPECT_EQ( TaskEnumToString( TaskEnum::SOLVE_FIELD ), "Solve" );

    auto& reg = TaskRegistry::Instance();
    reg.Register( "Solve", []() {
        return std::make_unique<SolveFieldContractTask>();
        } );

    auto byEnum = reg.Create( TaskEnum::SOLVE_FIELD );
    auto byName = reg.Create( "Solve" );
    ASSERT_NE( byEnum, nullptr );
    ASSERT_NE( byName, nullptr );
    EXPECT_TRUE( byEnum->NeedsSystemMap() );
    EXPECT_TRUE( byName->NeedsSystemMap() );
}

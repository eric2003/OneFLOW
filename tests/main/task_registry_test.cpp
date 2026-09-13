#include <gtest/gtest.h>
#include "SimuTask.h"
#include <algorithm>
#include <string>
#include <vector>

using namespace ONEFLOW;

namespace {

// Minimal mock task for registry tests.
class MockTask : public ISimuTask
{
public:
    explicit MockTask( bool needsMap = false )
        : needsMap_( needsMap )
    {
    }

    void Execute() override { executed_ = true; }
    bool NeedsSystemMap() const override { return needsMap_; }

    bool Executed() const { return executed_; }

private:
    bool needsMap_  = false;
    bool executed_  = false;
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

    task->Execute();
    // Execute ran without throw; mock has no public side-channel after move,
    // so just assert Create returned a valid object that can Execute.
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
    const std::string name = "Solve"; // may already exist if full binary linked;
                                      // for this test binary only mock registrations exist.

    // Ensure a known key is present for the enum path.
    reg.Register( "Solve", []() {
        return std::make_unique<MockTask>( true );
    } );

    auto task = reg.Create( TaskEnum::SOLVE_FIELD );
    ASSERT_NE( task, nullptr );
    EXPECT_TRUE( task->NeedsSystemMap() );
}

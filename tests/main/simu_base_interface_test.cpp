#include <gtest/gtest.h>
#include "SimuBase.h"
#include <memory>
#include <string>

using namespace ONEFLOW;

namespace {

// Stand-in for SimuImp: proves any SimuBase-derived type can be held and Run().
class FakeFullSimu : public SimuBase
{
public:
    void Run() override { ran_ = true; }
    bool Ran() const { return ran_; }

private:
    bool ran_ = false;
};

class FakeLightTest : public SimuBase
{
public:
    void Run() override { ran_ = true; }
    bool Ran() const { return ran_; }

private:
    bool ran_ = false;
};

} // namespace

TEST( SimuBaseInterfaceTest, FullAndLightShareSameInterface )
{
    std::unique_ptr<SimuBase> full  = std::make_unique<FakeFullSimu>();
    std::unique_ptr<SimuBase> light = std::make_unique<FakeLightTest>();

    ASSERT_NE( full, nullptr );
    ASSERT_NE( light, nullptr );

    full->Run();
    light->Run();

    // Downcast only for verification of side effect.
    EXPECT_TRUE( static_cast<FakeFullSimu*>( full.get() )->Ran() );
    EXPECT_TRUE( static_cast<FakeLightTest*>( light.get() )->Ran() );
}

TEST( SimuBaseInterfaceTest, TestRegistryCreatePolymorphic )
{
    auto& reg = TestRegistry::Instance();
    reg.Register( "fake_light", []() {
        return std::make_unique<FakeLightTest>();
    } );

    auto ptr = reg.Create( "fake_light" );
    ASSERT_NE( ptr, nullptr );
    ptr->Run();
    EXPECT_TRUE( static_cast<FakeLightTest*>( ptr.get() )->Ran() );
}

TEST( SimuBaseInterfaceTest, TestRegistryUnknownReturnsNull )
{
    auto ptr = TestRegistry::Instance().Create( "no_such_case_ever" );
    EXPECT_EQ( ptr, nullptr );
}

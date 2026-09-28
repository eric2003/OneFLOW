#include <gtest/gtest.h>
#include "WallDistPolicy.h"

using namespace ONEFLOW;

TEST( WallDistPolicy, Plateuns2dslau2_LaminarSkips )
{
    // vismodel=1, startStrategy=0, ireadwdst=0
    EXPECT_EQ(
        DecideFlowWallDistAction( 1, 0, 0 ),
        FlowWallDistAction::SkipAfterAlloc );
}

TEST( WallDistPolicy, EulerSkips )
{
    EXPECT_EQ(
        DecideFlowWallDistAction( 0, 0, 0 ),
        FlowWallDistAction::SkipAfterAlloc );
}

TEST( WallDistPolicy, TurbCreateWhenIreadZero )
{
    EXPECT_EQ(
        DecideFlowWallDistAction( 3, 0, 0 ),
        FlowWallDistAction::Create );
}

TEST( WallDistPolicy, TurbLoadWhenIreadNonZero )
{
    EXPECT_EQ(
        DecideFlowWallDistAction( 3, 0, 1 ),
        FlowWallDistAction::Load );
}

TEST( WallDistPolicy, RestartPrefersLoad )
{
    // startStrategy > 0 -> load even if ireadwdst == 0
    EXPECT_EQ(
        DecideFlowWallDistAction( 3, 1, 0 ),
        FlowWallDistAction::Load );
}
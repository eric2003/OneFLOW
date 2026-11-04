#include <gtest/gtest.h>
#include "TimeIntegralPolicy.h"
using namespace ONEFLOW;

TEST( TimeIntegralPolicy, RungeKuttaStageWeightMatchesOneOverStagePlusOne )
{
    EXPECT_DOUBLE_EQ( RungeKuttaStageWeight( 0 ), 1.0 );
    EXPECT_DOUBLE_EQ( RungeKuttaStageWeight( 1 ), 0.5 );
    EXPECT_DOUBLE_EQ( RungeKuttaStageWeight( 2 ), 1.0 / 3.0 );
    EXPECT_DOUBLE_EQ( RungeKuttaStageWeight( 3 ), 0.25 );
}

// Intermediate: stage weights on a stage-index grid.
TEST( TimeIntegralPolicy, RungeKuttaStageWeightCurveMatchesClosedForm )
{
    for ( int s = 0; s < 12; ++ s )
    {
        EXPECT_NEAR(
            RungeKuttaStageWeight( s ),
            1.0 / static_cast< Real >( s + 1 ),
            1.0e-12 );
    }
}

TEST( TimeIntegralPolicy, RungeKuttaStageWeightMonotoneDecreasing )
{
    for ( int s = 0; s < 8; ++ s )
    {
        EXPECT_GT( RungeKuttaStageWeight( s ), RungeKuttaStageWeight( s + 1 ) );
    }
}

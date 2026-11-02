#include <gtest/gtest.h>

#include "TurbSaPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( SaUpdateBound, UpperBoundFormula )
{
    // 10 * 2 / (1 * 0.5) = 40
    EXPECT_DOUBLE_EQ(
        SaWorkingVariableUpperBound( 10.0, 2.0, 1.0, 0.5 ), 40.0 );
}

TEST( SaUpdateBound, NonNegativeKeepsPositive )
{
    EXPECT_DOUBLE_EQ( SaNonNegativeWorkingVariable( 1.5, 0.1 ), 1.5 );
}

TEST( SaUpdateBound, NonNegativeReplacesNegative )
{
    EXPECT_DOUBLE_EQ( SaNonNegativeWorkingVariable( -2.0, 0.3 ), 0.3 );
}

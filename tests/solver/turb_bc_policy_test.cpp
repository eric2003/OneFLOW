#include <gtest/gtest.h>

#include "TurbBcPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TurbBcPolicy, SaWallImageNegates )
{
    EXPECT_DOUBLE_EQ( TurbWallSaImage( 3.0 ), -3.0 );
}

TEST( TurbBcPolicy, WallOmegaMenter )
{
    // 60 * 2 / (1 * 0.075 * 4 * 1) = 120 / 0.3 = 400
    EXPECT_DOUBLE_EQ(
        TurbWallOmegaMenter( 2.0, 1.0, 0.075, 4.0, 1.0 ), 400.0 );
}

TEST( TurbBcPolicy, KineticEnergyImage )
{
    EXPECT_DOUBLE_EQ( TurbWallKineticEnergyImage( 1.5 ), -1.5 );
}

TEST( TurbBcPolicy, FarfieldUsesInteriorWhenOutflow )
{
    EXPECT_TRUE( TurbFarfieldUsesInterior( 1.0 ) );
    EXPECT_TRUE( TurbFarfieldUsesInterior( 0.0 ) );
    EXPECT_FALSE( TurbFarfieldUsesInterior( -0.1 ) );
}

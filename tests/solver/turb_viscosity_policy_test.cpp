#include <gtest/gtest.h>
#include <limits>

#include "TurbViscosityPolicy.h"

using namespace ONEFLOW;

TEST( TurbViscosityPolicy, CapsEddyViscosityAtLaminarRatio )
{
    EXPECT_DOUBLE_EQ(
        TurbCappedEddyViscosity( 100.0, 2.0, 10.0 ), 20.0 );
}

TEST( TurbViscosityPolicy, PreservesEddyViscosityBelowLaminarRatio )
{
    EXPECT_DOUBLE_EQ(
        TurbCappedEddyViscosity( 5.0, 2.0, 10.0 ), 5.0 );
}
TEST( TurbViscosityPolicy, RejectsNegativeDensity )
{
    EXPECT_TRUE( TurbDensityIsInvalid( -1.0 ) );
}

TEST( TurbViscosityPolicy, RejectsNanDensity )
{
    EXPECT_TRUE( TurbDensityIsInvalid( std::numeric_limits< Real >::quiet_NaN() ) );
}

TEST( TurbViscosityPolicy, AcceptsZeroAndPositiveDensity )
{
    EXPECT_FALSE( TurbDensityIsInvalid( 0.0 ) );
    EXPECT_FALSE( TurbDensityIsInvalid( 1.2 ) );
}

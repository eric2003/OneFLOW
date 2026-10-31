#include <gtest/gtest.h>
#include <cmath>

#include "TurbSstPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TurbSstEddyViscosity, F2ArgTakesMax )
{
    EXPECT_DOUBLE_EQ( SstF2Arg( 1.0, 3.0 ), 3.0 );
}

TEST( TurbSstEddyViscosity, F2FromArgAtZeroIsZero )
{
    EXPECT_DOUBLE_EQ( SstF2FromArg( 0.0 ), 0.0 );
}

TEST( TurbSstEddyViscosity, F2FromArgLargeApproachesOne )
{
    EXPECT_NEAR( SstF2FromArg( 3.0 ), 1.0, 1.0e-6 );
}

TEST( TurbSstEddyViscosity, EddyViscosityLimitedByOmega )
{
    // vort*f2/a1 small => denominator = omega
    EXPECT_DOUBLE_EQ(
        SstEddyViscosity( 1.0, 2.0, 4.0, 0.0, 1.0, 0.31, 1.0 ), 0.5 );
}

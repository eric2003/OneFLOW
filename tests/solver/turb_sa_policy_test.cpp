#include <gtest/gtest.h>
#include <cmath>

#include "TurbSaPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TurbSaPolicy, SaChiIsProduct )
{
    EXPECT_NEAR( SaChi( 2.0, 0.5 ), 1.0 + SMALL, 1.0e-12 );
}

TEST( TurbSaPolicy, SaFv1AtLargeChiApproachesOne )
{
    const Real chi = 100.0;
    const Real cv13 = 7.1 * 7.1 * 7.1;
    EXPECT_NEAR( SaFv1( chi, cv13 ), 1.0, 1.0e-3 );
}

TEST( TurbSaPolicy, SaFv2AtZeroChiIsOne )
{
    EXPECT_DOUBLE_EQ( SaFv2( 0.0, 0.0 ), 1.0 );
}

TEST( TurbSaPolicy, SaModifiedStrainSimpleBranch )
{
    EXPECT_DOUBLE_EQ( SaModifiedStrainRate( 2.0, 1.0, 0.7, 0.9 ), 3.0 );
}

TEST( TurbSaPolicy, SaRRatioIsClamped )
{
    EXPECT_DOUBLE_EQ( SaRRatio( 100.0, 1.0 ), 10.0 );
}

TEST( TurbSaPolicy, SaGAndFwAtRZero )
{
    EXPECT_DOUBLE_EQ( SaG( 0.0, 0.3 ), 0.0 );
}

TEST( TurbSaPolicy, SaFt2DisabledIsZero )
{
    EXPECT_DOUBLE_EQ( SaFt2( 1.0, 1.2, 0.5, false ), 0.0 );
}

#include <gtest/gtest.h>

#include "TurbSstPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TurbSstPolicy, ProdkStrainOnly )
{
    EXPECT_DOUBLE_EQ( SstProdkStrainOnly( 2.0, 3.0, 0.5 ), 3.0 );
}

TEST( TurbSstPolicy, ProdkBoussinesqZeroDiv )
{
    EXPECT_DOUBLE_EQ(
        SstProdkBoussinesq( 2.0, 4.0, 0.0, 1.0, 1.0, 1.0 ), 8.0 );
}

TEST( TurbSstPolicy, DisskStandard )
{
    EXPECT_DOUBLE_EQ(
        SstDisskStandard( 1.0, 0.09, 1.0, 2.0, 3.0 ), 0.09 * 2.0 * 3.0 );
}

TEST( TurbSstPolicy, LimitProdkCaps )
{
    EXPECT_DOUBLE_EQ( SstLimitProdk( 10.0, 2.0, 3.0 ), 6.0 );
}

TEST( TurbSstPolicy, CrossDiffScalar )
{
    EXPECT_DOUBLE_EQ(
        SstCrossDiffScalar( 1.0, 0.5, 2.0, 1.0 ), 2.0 );
}

TEST( TurbSstPolicy, CdkwWithFloor )
{
    EXPECT_DOUBLE_EQ( SstCdkwWithFloor( 1.0e-30, 1.0e-20 ), 1.0e-20 );
}

TEST( TurbSstPolicy, BlendingArg1TakesMinOfMax )
{
    EXPECT_DOUBLE_EQ( SstBlendingArg1( 1.0, 3.0, 2.0 ), 2.0 );
    EXPECT_DOUBLE_EQ( SstBlendingArg1( 1.0, 0.5, 2.0 ), 1.0 );
}

TEST( TurbSstPolicy, F1FromArg1AtZeroIsZero )
{
    EXPECT_DOUBLE_EQ( SstF1FromArg1( 0.0 ), 0.0 );
}

TEST( TurbSstPolicy, F1FromArg1LargeApproachesOne )
{
    EXPECT_NEAR( SstF1FromArg1( 2.0 ), 1.0, 1.0e-6 );
}

TEST( TurbSstPolicy, GradientDotProduct )
{
    EXPECT_DOUBLE_EQ(
        SstGradientDotProduct( 1.0, 0.0, 0.0, 2.0, 0.0, 0.0 ), 2.0 );
}

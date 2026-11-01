#include <gtest/gtest.h>
#include <cmath>

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

TEST( TurbSstPolicy, BlendIsConvexCombination )
{
    // phi = F1*phi1 + (1-F1)*phi2
    EXPECT_DOUBLE_EQ( SstBlend( 0.25, 4.0, 8.0 ), 7.0 );
    EXPECT_DOUBLE_EQ( SstBlend( 1.0, 4.0, 8.0 ), 4.0 );
    EXPECT_DOUBLE_EQ( SstBlend( 0.0, 4.0, 8.0 ), 8.0 );
}

// Intermediate verification: F1(arg1) vs closed form tanh(arg1^4).
TEST( TurbSstPolicy, F1Arg1CurveMatchesClosedForm )
{
    const Real args[] = { 0.0, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0 };
    for ( Real a : args )
    {
        const Real expected = std::tanh( std::pow( a, 4 ) );
        EXPECT_NEAR( SstF1FromArg1( a ), expected, 1.0e-12 );
    }
}

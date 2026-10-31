#include <gtest/gtest.h>

#include "TurbSstPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TurbSstPolicy, ProdkStrainOnly )
{
    EXPECT_DOUBLE_EQ( SstProdkStrainOnly( 2.0, 3.0, 0.5 ), 3.0 );
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

TEST( TurbSstPolicy, BlendIsConvex )
{
    // F1*phi1 + (1-F1)*phi2 = 0.25*4 + 0.75*8 = 7
    EXPECT_DOUBLE_EQ( SstBlend( 0.25, 4.0, 8.0 ), 7.0 );
}

TEST( TurbSstPolicy, ProdwMenter )
{
    EXPECT_DOUBLE_EQ(
        SstProdwMenter( 0.5, 2.0, 1.0, 4.0, 1.0 ), 0.5 * 2.0 * 4.0 );
}

TEST( TurbSstPolicy, ProdwWilcox )
{
    EXPECT_DOUBLE_EQ(
        SstProdwWilcox( 2.0, 3.0, 1.0, 4.0 ), 2.0 * 3.0 * 4.0 );
}

TEST( TurbSstPolicy, Dissw )
{
    EXPECT_DOUBLE_EQ( SstDissw( 1.0, 0.07, 1.0, 2.0 ), 0.07 * 4.0 );
}

TEST( TurbSstPolicy, CrossDiffMenterZeroWhenF1IsOne )
{
    EXPECT_DOUBLE_EQ(
        SstCrossDiffMenter( 1.0, 1.0, 0.856, 2.0, 1.0 ), 0.0 );
}

TEST( TurbSstPolicy, CrossDiffWilcoxIgnoresNegative )
{
    EXPECT_DOUBLE_EQ(
        SstCrossDiffWilcox( 0.5, 1.0, -3.0, 2.0 ), 0.0 );
}

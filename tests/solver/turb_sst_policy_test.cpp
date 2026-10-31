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

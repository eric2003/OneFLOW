#include <gtest/gtest.h>

#include "TurbSstPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TurbSstResidual, SourceKIsProdMinusDiss )
{
    EXPECT_DOUBLE_EQ( SstSourceK( 5.0, 2.0 ), 3.0 );
}

TEST( TurbSstResidual, SourceWIncludesCrossDiff )
{
    EXPECT_DOUBLE_EQ( SstSourceW( 4.0, 1.0, 0.5 ), 3.5 );
}

TEST( TurbSstResidual, DiagonalK )
{
    EXPECT_DOUBLE_EQ( SstDiagonalK( 1.0, 0.09, 2.0 ), -0.36 );
}

TEST( TurbSstResidual, NegativeSourcePart )
{
    EXPECT_DOUBLE_EQ( SstNegativeSourcePart( 3.0 ), 0.0 );
    EXPECT_DOUBLE_EQ( SstNegativeSourcePart( -4.0 ), -4.0 );
}

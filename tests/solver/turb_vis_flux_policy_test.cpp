#include <gtest/gtest.h>

#include "TurbVisFluxPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TurbVisFluxPolicy, ArithmeticMean )
{
    EXPECT_DOUBLE_EQ( ArithmeticMean( 2.0, 4.0 ), 3.0 );
}

TEST( TurbVisFluxPolicy, SaFaceEffectiveDiffusivity )
{
    EXPECT_DOUBLE_EQ( SaFaceEffectiveDiffusivity( 1.0, 1.0, 0.5 ), 1.0 );
}

TEST( TurbVisFluxPolicy, SstFaceEffectiveDiffusivity )
{
    EXPECT_DOUBLE_EQ( SstFaceEffectiveDiffusivity( 2.0, 3.0, 0.5 ), 3.5 );
}

TEST( TurbVisFluxPolicy, TurbScalarDiffusiveFluxSign )
{
    EXPECT_DOUBLE_EQ(
        TurbScalarDiffusiveFlux( 0.5, 2.0, 3.0, 4.0 ), -12.0 );
}

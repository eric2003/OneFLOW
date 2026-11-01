#include <gtest/gtest.h>

#include "NsLusgsPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( LusgsViscousSpectrum, ContributionIsMinusRadiusTimesDq )
{
    EXPECT_DOUBLE_EQ(
        LusgsViscousSpectrumContribution( 2.0, 3.0 ), -6.0 );
}

TEST( LusgsViscousSpectrum, FaceArithmeticMean )
{
    EXPECT_DOUBLE_EQ( LusgsFaceArithmeticMean( 2.0, 4.0 ), 3.0 );
}

TEST( LusgsViscousSpectrum, RadiusPerVolume )
{
    EXPECT_DOUBLE_EQ( LusgsViscousRadiusPerVolume( 10.0, 2.0 ), 5.0 );
}

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

// Intermediate: d(rhs) = -visrad * dq on a dq grid (visrad=2).
TEST( LusgsViscousSpectrum, ContributionCurveMatchesClosedForm )
{
    const Real visrad = 2.0;
    const Real dqs[] = { -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0 };
    for ( Real dq : dqs )
    {
        EXPECT_NEAR(
            LusgsViscousSpectrumContribution( visrad, dq ),
            -visrad * dq,
            1.0e-12 );
    }
}

// Intermediate: face mean and per-volume radius samples.
TEST( LusgsViscousSpectrum, MeanAndRadiusClosedForm )
{
    EXPECT_DOUBLE_EQ( LusgsFaceArithmeticMean( 1.0, 3.0 ), 2.0 );
    EXPECT_DOUBLE_EQ( LusgsViscousRadiusPerVolume( 8.0, 4.0 ), 2.0 );
}

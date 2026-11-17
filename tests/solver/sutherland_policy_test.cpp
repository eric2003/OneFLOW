#include <gtest/gtest.h>
#include <cmath>

#include "Sutherland.h"

using namespace ONEFLOW;

TEST( SutherlandPolicy, NondimViscosityAtUnityIsOne )
{
    // T=1, ts=0.4 -> (1)^{3/2} * 1.4 / 1.4 = 1
    EXPECT_NEAR( SutherlandViscosityNondim( 1.0, 0.4 ), 1.0, 1.0e-12 );
}

TEST( SutherlandPolicy, NondimViscosityMatchesFormula )
{
    const Real t = 2.0;
    const Real ts = 0.5;
    const Real expected = std::sqrt( t * t * t ) * ( 1.0 + ts ) / ( t + ts );
    EXPECT_NEAR( SutherlandViscosityNondim( t, ts ), expected, 1.0e-12 );
}

TEST( SutherlandPolicy, DimViscosityMatchesClassicSutherland )
{
    const Real t0 = 273.16;
    const Real ts = 110.4;
    const Real mu0 = 1.715e-5;
    const Real t = 300.0;
    const Real expected =
        ( ( t0 + ts ) / ( t + ts ) ) * std::pow( t / t0, 1.5 ) * mu0;
    EXPECT_NEAR( SutherlandViscosityDim( t, t0, ts, mu0 ), expected, 1.0e-18 );
}

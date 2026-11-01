#include <gtest/gtest.h>

#include "TimeStepPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TimeStepPolicy, InviscidFaceSpectralRadiusMatchesFormula )
{
    EXPECT_DOUBLE_EQ(
        InviscidFaceSpectralRadius( -3.0, 1.0, 4.0 ),
        half * 4.0 * ( 3.0 + 1.0 ) );
}

TEST( TimeStepPolicy, ViscousFaceSpectralRadiusModel3IsArea2OverReRho )
{
    EXPECT_DOUBLE_EQ(
        ViscousFaceSpectralRadiusModel3( 2.0, 100.0, 3.0 ),
        9.0 / ( 100.0 * 2.0 ) );
}

TEST( TimeStepPolicy, CellTimeStepFromSpectrumIsCflVolOverLambda )
{
    EXPECT_DOUBLE_EQ( CellTimeStepFromSpectrum( 0.5, 4.0, 2.0 ), 1.0 );
}

TEST( TimeStepPolicy, CombineInvAndVisTimeStepIsHarmonic )
{
    EXPECT_DOUBLE_EQ( CombineInvAndVisTimeStep( 2.0, 2.0 ), 1.0 );
}

// Intermediate: lambda = (A/2)(|vn|+c) on a |vn| grid (A=2,c=1).
TEST( TimeStepPolicy, InviscidSpectrumCurveMatchesClosedForm )
{
    const Real faceArea = 2.0;
    const Real sound = 1.0;
    const Real absVns[] = { 0.0, 0.5, 1.0, 2.0, 5.0, 10.0 };
    for ( Real absVn : absVns )
    {
        // vnRelMid can be +/- absVn
        const Real expected = half * faceArea * ( absVn + sound );
        EXPECT_NEAR(
            InviscidFaceSpectralRadius( absVn, sound, faceArea ),
            expected,
            1.0e-12 );
        EXPECT_NEAR(
            InviscidFaceSpectralRadius( -absVn, sound, faceArea ),
            expected,
            1.0e-12 );
    }
}

// Intermediate: dt = CFL * Vol / lambda
TEST( TimeStepPolicy, CellTimeStepCurveMatchesClosedForm )
{
    const Real cfl = 1.0;
    const Real vol = 2.0;
    const Real lams[] = { 0.5, 1.0, 2.0, 4.0 };
    for ( Real lam : lams )
    {
        EXPECT_NEAR(
            CellTimeStepFromSpectrum( cfl, vol, lam ),
            cfl * vol / lam,
            1.0e-12 );
    }
}

// Intermediate: harmonic combine dt_inv * dt_vis / (dt_inv + dt_vis)
TEST( TimeStepPolicy, CombineHarmonicMoreSamples )
{
    EXPECT_NEAR( CombineInvAndVisTimeStep( 1.0, 3.0 ), 0.75, 1.0e-12 );
    EXPECT_NEAR( CombineInvAndVisTimeStep( 4.0, 4.0 ), 2.0, 1.0e-12 );
}

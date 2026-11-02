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

// Intermediate: ViscousFaceSpectralRadiusModel1 on a mut grid.
// Ref: rho=1, mu=1, invPr_l=1/0.72, invPr_t=1/0.9, gamma=1.4, Re=100, A=1.
TEST( TimeStepPolicy, ViscousFaceSpectralRadiusModel1CurveMatchesClosedForm )
{
    const Real density = 1.0;
    const Real mu = 1.0;
    const Real invPrL = 1.0 / 0.72;
    const Real invPrT = 1.0 / 0.9;
    const Real gama = 1.4;
    const Real reynolds = 100.0;
    const Real area = 1.0;
    const Real muts[] = { 0.0, 0.5, 1.0, 2.0, 5.0, 10.0 };
    for ( Real mut : muts )
    {
        const Real c1 = ( 4.0 / 3.0 ) * ( mu + mut );
        const Real c2 = gama * ( mu * invPrL + mut * invPrT );
        const Real c3 = two * MAX( c1, c2 ) / ( reynolds * density );
        const Real expected = SQR( area ) * c3;
        EXPECT_NEAR(
            ViscousFaceSpectralRadiusModel1(
                density, mu, mut, invPrL, invPrT, gama, reynolds, area ),
            expected,
            1.0e-12 );
    }
}

TEST( TimeStepPolicy, ViscousFaceSpectralRadiusModel1ScalesWithAreaSquared )
{
    const Real base = ViscousFaceSpectralRadiusModel1(
        1.0, 1.0, 1.0, 1.0 / 0.72, 1.0 / 0.9, 1.4, 100.0, 1.0 );
    const Real twice = ViscousFaceSpectralRadiusModel1(
        1.0, 1.0, 1.0, 1.0 / 0.72, 1.0 / 0.9, 1.4, 100.0, 2.0 );
    EXPECT_NEAR( twice, 4.0 * base, 1.0e-12 );
}

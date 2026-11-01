#include <gtest/gtest.h>
#include <cmath>

#include "TurbSstPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TurbSstEddyViscosity, StrainRateMagnitudeMatchesSymmetricTensorNorm )
{
    // A diagonal strain with S11=1 and S22=-1 has |S|=sqrt(2*(1+1))=2.
    EXPECT_DOUBLE_EQ(
        SstStrainRateMagnitude(
            1.0, 0.0, 0.0,
            0.0, -1.0, 0.0,
            0.0, 0.0, 0.0 ),
        2.0 );
}

TEST( TurbSstEddyViscosity, RateMagnitudesDistinguishShearFromRotation )
{
    // du/dy=3 and dv/dx=1 give S12=2, |S|=4 and |Omega|=2.
    const Real strainRate = SstStrainRateMagnitude(
        0.0, 3.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 0.0, 0.0 );
    const Real vorticity = SstVorticityMagnitude(
        0.0, 3.0, 0.0,
        1.0, 0.0, 0.0,
        0.0, 0.0, 0.0 );
    EXPECT_DOUBLE_EQ( strainRate, 4.0 );
    EXPECT_DOUBLE_EQ( vorticity, 2.0 );
}

TEST( TurbSstEddyViscosity, ProductionOptionSelectsStrainRate )
{
    EXPECT_DOUBLE_EQ(
        SstLimitingRateMagnitude(
            1,
            0.0, 3.0, 0.0,
            1.0, 0.0, 0.0,
            0.0, 0.0, 0.0 ),
        4.0 );
}

TEST( TurbSstEddyViscosity, OtherProductionOptionsSelectVorticity )
{
    const Real expectedVorticity = 2.0;
    const int productionOptions[] = { 0, 2, -1 };
    for ( int productionOption : productionOptions )
    {
        EXPECT_DOUBLE_EQ(
            SstLimitingRateMagnitude(
                productionOption,
                0.0, 3.0, 0.0,
                1.0, 0.0, 0.0,
                0.0, 0.0, 0.0 ),
            expectedVorticity );
    }
}

TEST( TurbSstEddyViscosity, F2ArgFromTermsTakesMax )
{
    EXPECT_DOUBLE_EQ( SstF2ArgFromTerms( 1.0, 3.0 ), 3.0 );
}

TEST( TurbSstEddyViscosity, F2ArgumentUsesLengthScaleTermWhenLarger )
{
    // 2*sqrt(k)/(betaStar*omega*d) = 200/9; viscous term = 0.0025.
    EXPECT_NEAR(
        SstF2Argument( 4.0, 2.0, 1.0, 0.09, 0.001, 1.0, 100.0 ),
        200.0 / 9.0,
        1.0e-12 );
}

TEST( TurbSstEddyViscosity, F2ArgumentUsesViscousTermWhenLarger )
{
    // Length-scale term = 1/45; viscous term = 500.
    EXPECT_DOUBLE_EQ(
        SstF2Argument( 1.0e-8, 1.0, 0.1, 0.09, 0.01, 1.0, 1.0 ),
        500.0 );
}

TEST( TurbSstEddyViscosity, F2FromArgAtZeroIsZero )
{
    EXPECT_DOUBLE_EQ( SstF2FromArg( 0.0 ), 0.0 );
}

TEST( TurbSstEddyViscosity, F2FromArgLargeApproachesOne )
{
    EXPECT_NEAR( SstF2FromArg( 3.0 ), 1.0, 1.0e-6 );
}

TEST( TurbSstEddyViscosity, EddyViscosityLimitedByOmega )
{
    // vort*f2/a1 small => denominator = omega
    EXPECT_DOUBLE_EQ(
        SstEddyViscosity( 1.0, 2.0, 4.0, 0.0, 1.0, 0.31, 1.0 ), 0.5 );
}

TEST( TurbSstEddyViscosity, EddyViscosityLimitedByStrainOrVorticity )
{
    // The selected rate magnitude activates the Bradshaw limiter:
    // R*F2/a1 = 8*0.5/0.25 = 16, greater than omega = 4.
    // mu_t = rho*k/denominator*Re = 2*3/16*10 = 3.75.
    EXPECT_DOUBLE_EQ(
        SstEddyViscosity( 2.0, 3.0, 4.0, 8.0, 0.5, 0.25, 10.0 ),
        3.75 );
}

TEST( TurbSstEddyViscosity, CompleteCellClosureReportsEveryStage )
{
    SstEddyViscosityInput input;
    input.kineticEnergy = 3.0;
    input.omega = 4.0;
    input.wallDistance = 1.0;
    input.betaStar = 0.09;
    input.laminarViscosity = 0.001;
    input.density = 1.0;
    input.reynolds = 10.0;
    input.a1 = 0.31;
    input.maxViscosityRatio = 10000.0;
    input.productionOption = 1;
    input.dudy = 3.0;
    input.dvdx = 1.0;

    const SstEddyViscosityResult result = SstComputeEddyViscosity( input );
    const Real expectedArg2 = 2.0 * std::sqrt( 3.0 ) / ( 0.09 * 4.0 );
    const Real expectedF2 = std::tanh( expectedArg2 * expectedArg2 );
    const Real expectedStrainRate = 4.0;
    const Real expectedRawViscosity =
        3.0 * 10.0 / std::max( 4.0, expectedStrainRate * expectedF2 / 0.31 );

    EXPECT_NEAR( result.f2Argument, expectedArg2, 1.0e-12 );
    EXPECT_NEAR( result.f2, expectedF2, 1.0e-12 );
    EXPECT_DOUBLE_EQ( result.limitingRateMagnitude, expectedStrainRate );
    EXPECT_NEAR( result.uncappedEddyViscosity, expectedRawViscosity, 1.0e-12 );
    EXPECT_NEAR( result.eddyViscosity, expectedRawViscosity, 1.0e-12 );
}

TEST( TurbSstEddyViscosity, CompleteCellClosureUsesVorticityAndSharedCap )
{
    SstEddyViscosityInput input;
    input.kineticEnergy = 3.0;
    input.omega = 4.0;
    input.wallDistance = 1.0;
    input.betaStar = 0.09;
    input.laminarViscosity = 0.001;
    input.density = 1.0;
    input.reynolds = 10.0;
    input.a1 = 0.31;
    input.maxViscosityRatio = 1.0;
    input.productionOption = 0;
    input.dudy = 3.0;
    input.dvdx = 1.0;

    const SstEddyViscosityResult result = SstComputeEddyViscosity( input );
    const Real expectedVorticity = 2.0;
    const Real expectedArg2 = 2.0 * std::sqrt( 3.0 ) / ( 0.09 * 4.0 );
    const Real expectedF2 = std::tanh( expectedArg2 * expectedArg2 );
    const Real expectedRawViscosity =
        3.0 * 10.0 / std::max( 4.0, expectedVorticity * expectedF2 / 0.31 );

    EXPECT_DOUBLE_EQ( result.limitingRateMagnitude, expectedVorticity );
    EXPECT_NEAR( result.uncappedEddyViscosity, expectedRawViscosity, 1.0e-12 );
    EXPECT_DOUBLE_EQ( result.eddyViscosity, 0.001 );
}

// Intermediate verification: F2(arg2) vs closed form tanh(arg2^2).
// Proves pure SstF2FromArg; does not prove full SST eddy viscosity on a mesh.
TEST( TurbSstEddyViscosity, F2Arg2CurveMatchesClosedForm )
{
    const Real args[] = { 0.0, 0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 3.0 };
    for ( Real a : args )
    {
        const Real expected = std::tanh( a * a );
        EXPECT_NEAR( SstF2FromArg( a ), expected, 1.0e-12 );
    }
}
 
TEST( TurbSstEddyViscosity, StrainAndVorticityMagnitudesAreRotationInvariant )
{
    const Real gradient[3][3] = {
        { 1.0, 2.0, -0.5 },
        { -3.0, 0.25, 1.5 },
        { 0.75, -1.25, -0.5 }
    };
    const Real rotation[3][3] = {
        { 0.0, -1.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 0.0, 0.0, 1.0 }
    };
    Real rotated[3][3] = {};
    Real intermediate[3][3] = {};

    // Change basis with G' = Q^T G Q; both tensor norms must be invariant.
    for ( int i = 0; i < 3; ++i )
    {
        for ( int j = 0; j < 3; ++j )
        {
            for ( int k = 0; k < 3; ++k )
            {
                intermediate[i][j] += gradient[i][k] * rotation[k][j];
            }
        }
    }
    for ( int i = 0; i < 3; ++i )
    {
        for ( int j = 0; j < 3; ++j )
        {
            for ( int k = 0; k < 3; ++k )
            {
                rotated[i][j] += rotation[k][i] * intermediate[k][j];
            }
        }
    }

    const Real strain = SstStrainRateMagnitude(
        gradient[0][0], gradient[0][1], gradient[0][2],
        gradient[1][0], gradient[1][1], gradient[1][2],
        gradient[2][0], gradient[2][1], gradient[2][2] );
    const Real rotatedStrain = SstStrainRateMagnitude(
        rotated[0][0], rotated[0][1], rotated[0][2],
        rotated[1][0], rotated[1][1], rotated[1][2],
        rotated[2][0], rotated[2][1], rotated[2][2] );
    const Real vorticity = SstVorticityMagnitude(
        gradient[0][0], gradient[0][1], gradient[0][2],
        gradient[1][0], gradient[1][1], gradient[1][2],
        gradient[2][0], gradient[2][1], gradient[2][2] );
    const Real rotatedVorticity = SstVorticityMagnitude(
        rotated[0][0], rotated[0][1], rotated[0][2],
        rotated[1][0], rotated[1][1], rotated[1][2],
        rotated[2][0], rotated[2][1], rotated[2][2] );

    EXPECT_NEAR( rotatedStrain, strain, 1.0e-12 );
    EXPECT_NEAR( rotatedVorticity, vorticity, 1.0e-12 );
}

TEST( TurbSstEddyViscosity, SolidBodyRotationHasNoStrain )
{
    // A rigid rotation has an antisymmetric velocity gradient.
    const Real strain = SstStrainRateMagnitude(
        0.0, -2.0, 0.0,
        2.0, 0.0, 0.0,
        0.0, 0.0, 0.0 );
    const Real vorticity = SstVorticityMagnitude(
        0.0, -2.0, 0.0,
        2.0, 0.0, 0.0,
        0.0, 0.0, 0.0 );

    EXPECT_DOUBLE_EQ( strain, 0.0 );
    EXPECT_DOUBLE_EQ( vorticity, 4.0 );
}


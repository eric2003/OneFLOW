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

TEST( TurbSstEddyViscosity, F2ArgTakesMax )
{
    EXPECT_DOUBLE_EQ( SstF2Arg( 1.0, 3.0 ), 3.0 );
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

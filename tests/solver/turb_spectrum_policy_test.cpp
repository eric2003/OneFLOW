#include <gtest/gtest.h>

#include "TurbSpectrumPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TurbSpectrumPolicy, InviscidFaceSpectralRadius )
{
    EXPECT_DOUBLE_EQ( TurbInviscidFaceSpectralRadius( 4.0, 2.0 ), 4.0 );
}

TEST( TurbSpectrumPolicy, InviscidFaceMatrixContribSigns )
{
    Real ml, mr;
    TurbInviscidFaceMatrixContrib( 2.0, 1.0, -1.0, 1.0, ml, mr );
    EXPECT_DOUBLE_EQ( ml, 2.0 * half * ( -1.0 - 1.0 ) );
    EXPECT_DOUBLE_EQ( mr, 2.0 * half * ( -1.0 - 1.0 ) );
}

TEST( TurbSpectrumPolicy, SaEffectiveKinematicViscosity )
{
    EXPECT_DOUBLE_EQ( SaEffectiveKinematicViscosity( 2.0, 1.0, 0.5 ), 2.5 );
}

TEST( TurbSpectrumPolicy, DiffusiveFaceSpectralRadius )
{
    EXPECT_DOUBLE_EQ(
        TurbDiffusiveFaceSpectralRadius( 0.5, 2.0, 3.0, 4.0, 0.25 ),
        0.5 * 0.25 * ( 2.0 * 3.0 * 4.0 ) );
}

// Intermediate: lambda_inv = (A/2) * |vn| on a |vn| grid (A=2).
TEST( TurbSpectrumPolicy, InviscidSpectrumCurveMatchesClosedForm )
{
    const Real faceArea = 2.0;
    const Real absVns[] = { 0.0, 0.5, 1.0, 2.0, 5.0, 10.0 };
    for ( Real absVn : absVns )
    {
        const Real expected = faceArea * half * absVn;
        EXPECT_NEAR(
            TurbInviscidFaceSpectralRadius( faceArea, absVn ),
            expected,
            1.0e-12 );
    }
}

// Intermediate: nu_et = mu_l/rho + nu_tilde
TEST( TurbSpectrumPolicy, SaEffectiveNuMatchesClosedForm )
{
    EXPECT_DOUBLE_EQ( SaEffectiveKinematicViscosity( 1.0, 2.0, 0.25 ), 0.75 );
    EXPECT_DOUBLE_EQ( SaEffectiveKinematicViscosity( 0.0, 1.0, 3.0 ), 3.0 );
}

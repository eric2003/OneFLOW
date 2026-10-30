#include <gtest/gtest.h>
#include <cmath>

#include "NsBcPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( BcPolicy, ReflectVelocityCancelsNormalComponent )
{
    Real ur, vr, wr;
    ReflectVelocityRelativeToMovingFace( 3.0, 4.0, 0.0, 1.0, 0.0, 0.0, 0.0, ur, vr, wr );
    EXPECT_DOUBLE_EQ( ur, -3.0 );
    EXPECT_DOUBLE_EQ( vr, 4.0 );
}

TEST( BcPolicy, ImageVelocityAboutStationaryWallNegates )
{
    Real ui, vi, wi;
    ImageVelocityAboutWall( 2.0, -1.0, 0.5, 0.0, 0.0, 0.0, ui, vi, wi );
    EXPECT_DOUBLE_EQ( ui, -2.0 );
    EXPECT_DOUBLE_EQ( vi, 1.0 );
}

TEST( BcPolicy, ArithmeticAverageIsMidpoint )
{
    EXPECT_DOUBLE_EQ( ArithmeticAverage( 1.0, 3.0 ), 2.0 );
}

TEST( BcPolicy, SubsonicAcousticRiemannCombineMatchesInvariants )
{
    Real vnb, cb;
    // vn_in=0, c_in=1, vn_ref=0, c_ref=1, g-1=0.4
    // R+ = 0 + 2/0.4 = 5, R- = 0 - 5 = -5, vnb=0, cb = 0.25*10*0.4 = 1
    SubsonicAcousticRiemannCombine( 0.0, 1.0, 0.0, 1.0, 0.4, vnb, cb );
    EXPECT_NEAR( vnb, 0.0, 1.0e-12 );
    EXPECT_NEAR( cb, 1.0, 1.0e-12 );
}

TEST( BcPolicy, TangentialAndReconstructAreInverse )
{
    Real vtx, vty, vtz;
    TangentialVelocityFromNormalProjection( 3.0, 4.0, 0.0, 1.0, 0.0, 0.0, 3.0, vtx, vty, vtz );
    EXPECT_DOUBLE_EQ( vtx, 0.0 );
    EXPECT_DOUBLE_EQ( vty, 4.0 );
    Real u, v, w;
    ReconstructVelocityFromTangentialAndNormal( vtx, vty, vtz, 1.0, 0.0, 0.0, 3.0, u, v, w );
    EXPECT_DOUBLE_EQ( u, 3.0 );
    EXPECT_DOUBLE_EQ( v, 4.0 );
}

TEST( BcPolicy, DensityAndPressureFromIsentropicSound )
{
    // gamma=1.4, entropy=1, c=sqrt(1.4) -> rho = 1, p = 1
    const Real gama = 1.4;
    const Real c = std::sqrt( gama );
    const Real rho = DensityFromIsentropicSoundSpeed( c, 1.0, gama, gama - 1.0 );
    EXPECT_NEAR( rho, 1.0, 1.0e-12 );
    EXPECT_NEAR( PressureFromSoundSpeedDensity( c, rho, gama ), 1.0, 1.0e-12 );
}

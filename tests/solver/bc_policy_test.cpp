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
    EXPECT_DOUBLE_EQ( wr, 0.0 );
}

TEST( BcPolicy, ImageVelocityAboutStationaryWallNegates )
{
    Real ui, vi, wi;
    ImageVelocityAboutWall( 2.0, -1.0, 0.5, 0.0, 0.0, 0.0, ui, vi, wi );
    EXPECT_DOUBLE_EQ( ui, -2.0 );
    EXPECT_DOUBLE_EQ( vi, 1.0 );
    EXPECT_DOUBLE_EQ( wi, -0.5 );
}

TEST( BcPolicy, ImageVelocityAboutMovingWall )
{
    Real ui, vi, wi;
    ImageVelocityAboutWall( 1.0, 0.0, 0.0, 3.0, 0.0, 0.0, ui, vi, wi );
    EXPECT_DOUBLE_EQ( ui, 5.0 );
}

TEST( BcPolicy, ArithmeticAverageIsMidpoint )
{
    EXPECT_DOUBLE_EQ( ArithmeticAverage( 1.0, 3.0 ), 2.0 );
}

TEST( BcPolicy, SubsonicAcousticRiemannCombineMatchesInvariants )
{
    Real vnb, cb;
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
    const Real gama = 1.4;
    const Real c = std::sqrt( gama );
    const Real rho = DensityFromIsentropicSoundSpeed( c, 1.0, gama, gama - 1.0 );
    EXPECT_NEAR( rho, 1.0, 1.0e-12 );
    EXPECT_NEAR( PressureFromSoundSpeedDensity( c, rho, gama ), 1.0, 1.0e-12 );
}

TEST( BcPolicy, IdealGasDensityFromPressureTemperature )
{
    EXPECT_DOUBLE_EQ(
        IdealGasDensityFromPressureTemperature( 2.0, 4.0, 0.5, 1.0 ), 1.0 );
}

TEST( BcPolicy, MirrorScalarWithFloorClamps )
{
    EXPECT_DOUBLE_EQ( MirrorScalarWithFloor( 10.0, 1.0, 0.5 ), 19.0 );
    EXPECT_DOUBLE_EQ( MirrorScalarWithFloor( 1.0, 10.0, 0.5 ), 0.5 );
}

TEST( BcPolicy, MirrorDensityAboutWallFloorsAtHalfWall )
{
    EXPECT_DOUBLE_EQ( MirrorDensityAboutWall( 2.0, 0.1 ), 3.9 );
    EXPECT_DOUBLE_EQ( MirrorDensityAboutWall( 2.0, 10.0 ), 1.0 );
}

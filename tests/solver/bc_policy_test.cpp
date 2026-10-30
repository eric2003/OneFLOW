#include <gtest/gtest.h>

#include "NsBcPolicy.h"

using namespace ONEFLOW;

TEST( BcPolicy, ReflectVelocityCancelsNormalComponent )
{
    Real ur, vr, wr;
    // n = (1,0,0), vfn=0, v=(3,4,0) -> vn=3, v'=( -3,4,0 )
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

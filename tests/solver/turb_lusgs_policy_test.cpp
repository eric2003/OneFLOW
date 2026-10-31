#include <gtest/gtest.h>

#include "TurbLusgsPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( TurbLusgsPolicy, MatrixTimesDq )
{
    EXPECT_DOUBLE_EQ( TurbLusgsMatrixTimesDq( 2.0, 3.0 ), 6.0 );
}

TEST( TurbLusgsPolicy, SolveByRadius )
{
    EXPECT_DOUBLE_EQ( TurbLusgsSolveByRadius( 10.0, 2.0 ), 5.0 );
}

TEST( TurbLusgsPolicy, Difference )
{
    EXPECT_DOUBLE_EQ( TurbLusgsDifference( 7.0, 3.0 ), 4.0 );
}

TEST( TurbLusgsPolicy, SweepNormContribution )
{
    EXPECT_DOUBLE_EQ( TurbLusgsSweepNormContribution( 5.0, 2.0 ), 9.0 );
}

#include "PointMachine.h"

#include <gtest/gtest.h>

namespace ONEFLOW
{

TEST( PointMachineLifecycleTest, OwnsPointsAndExposesNonOwningViews )
{
    PointMachine pointMachine;

    pointMachine.AddPoint( 1.0, 2.0, 3.0, 1 );
    pointMachine.AddPoint( 4.0, 5.0, 6.0, 2 );

    EXPECT_EQ( pointMachine.GetNPoint(), 2 );

    PointType & first = pointMachine.GetPoint( 1 );
    EXPECT_DOUBLE_EQ( first.x, 1.0 );
    EXPECT_DOUBLE_EQ( first.y, 2.0 );
    EXPECT_DOUBLE_EQ( first.z, 3.0 );

    const PointMachine & constPointMachine = pointMachine;
    const PointType & second = constPointMachine.GetPoint( 2 );

    EXPECT_DOUBLE_EQ( second.x, 4.0 );
    EXPECT_DOUBLE_EQ( second.y, 5.0 );
    EXPECT_DOUBLE_EQ( second.z, 6.0 );
}

TEST( PointMachineLifecycleTest, ResetReleasesAllOwnedPoints )
{
    PointMachine pointMachine;

    pointMachine.AddPoint( 1.0, 2.0, 3.0, 1 );
    pointMachine.Reset();

    EXPECT_EQ( pointMachine.GetNPoint(), 0 );
}

} // namespace ONEFLOW

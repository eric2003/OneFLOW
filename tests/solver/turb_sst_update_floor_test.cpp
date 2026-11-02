#include <gtest/gtest.h>

#include "TurbSstPolicy.h"
#include "HXMath.h"

using namespace ONEFLOW;

TEST( SstUpdateFloor, FloorKineticEnergy )
{
    EXPECT_DOUBLE_EQ( SstFloorKineticEnergy( 0.1, 1.0 ), 1.0 );
    EXPECT_DOUBLE_EQ( SstFloorKineticEnergy( 2.0, 1.0 ), 2.0 );
}

TEST( SstUpdateFloor, FloorOmega )
{
    EXPECT_DOUBLE_EQ( SstFloorOmega( 0.01, 0.05 ), 0.05 );
}

TEST( SstUpdateFloor, DensityScaledUpdate )
{
    EXPECT_DOUBLE_EQ( SstDensityScaledUpdate( 1.0, 4.0, 2.0 ), 3.0 );
}

TEST( SstUpdateFloor, TransitionCaps )
{
    EXPECT_DOUBLE_EQ( SstCapIntermittency( 1.5, 1.0 ), 1.0 );
    EXPECT_DOUBLE_EQ( SstFloorReTheta( 10.0, 20.0 ), 20.0 );
}

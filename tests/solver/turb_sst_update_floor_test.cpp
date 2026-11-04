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

// Intermediate: SstFloorOmega on an omega grid.
// Ref: floor=0.05; floored = max(omega, 0.05).
TEST( SstUpdateFloor, SstFloorOmegaCurveMatchesClosedForm )
{
    const Real floorValue = 0.05;
    const Real omegas[] = { 0.0, 0.01, 0.05, 0.1, 1.0 };
    for ( Real omega : omegas )
    {
        const Real expected = ( omega > floorValue ) ? omega : floorValue;
        EXPECT_NEAR( SstFloorOmega( omega, floorValue ), expected, 1.0e-12 );
    }
}

TEST( SstUpdateFloor, SstFloorOmegaNeverBelowFloor )
{
    EXPECT_GE( SstFloorOmega( -1.0, 0.05 ), 0.05 );
    EXPECT_GE( SstFloorOmega( 0.0, 0.05 ), 0.05 );
}

// Intermediate: SstDensityScaledUpdate on a deltaQ grid.
// Ref: q0=1, density=2; q = 1 + deltaQ/2.
TEST( SstUpdateFloor, SstDensityScaledUpdateCurveMatchesClosedForm )
{
    const Real q0 = 1.0;
    const Real density = 2.0;
    const Real dqs[] = { -4.0, -1.0, 0.0, 2.0, 4.0 };
    for ( Real dq : dqs )
    {
        EXPECT_NEAR(
            SstDensityScaledUpdate( q0, dq, density ),
            q0 + dq / density,
            1.0e-12 );
    }
}

TEST( SstUpdateFloor, SstDensityScaledUpdateInverselyScalesWithDensity )
{
    EXPECT_NEAR(
        SstDensityScaledUpdate( 0.0, 4.0, 1.0 ),
        2.0 * SstDensityScaledUpdate( 0.0, 4.0, 2.0 ),
        1.0e-12 );
}

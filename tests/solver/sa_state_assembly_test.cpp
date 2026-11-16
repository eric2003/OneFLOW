/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

// The field-to-view assembly that runs inside the solver, tested without a grid.
//
// Every device contract test in this repository builds its view by hand, so the
// kernels are verified while the code that fills the view in the solver is not.
// That gap was real: the SA velocity-gradient slot order was written transposed
// in UTurbSrcFlux::CalcSrcFlux1Equ - [dudx,dvdx,dwdx,...] instead of
// [dudx,dudy,dudz,...] - which compiles, passes every contract test, and
// produces wrong residuals silently.
//
// This test pins the mapping by giving each of the nine planes a distinct
// sentinel value and checking the exact slot each one lands in.

#include "SaStateAssembly.h"

#include <gtest/gtest.h>

#include <array>

using namespace ONEFLOW;

namespace
{

// Distinct, non-interchangeable values.  Using 1..9 rather than something
// symmetric matters: a transposed mapping would map 2->3 and 3->2, which any
// symmetric or repeated sentinel would hide.
struct GradientFixture
{
    std::array< std::array< Real, 3 >, 3 > storage = {};

    GradientFixture()
    {
        Real next = 1.0;
        for ( int component = 0; component < 3; ++ component )
        {
            for ( int direction = 0; direction < 3; ++ direction )
            {
                storage[ component ][ direction ] = next;
                next += 1.0;
            }
        }
    }

    SaMeanFlowGradientFields Fields()
    {
        SaMeanFlowGradientFields fields;
        for ( int component = 0; component < 3; ++ component )
        {
            for ( int direction = 0; direction < 3; ++ direction )
            {
                fields.gradient[ component ][ direction ] =
                    & storage[ component ][ direction ];
            }
        }
        return fields;
    }
};

}

// Slot order is dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz.
TEST( SaStateAssembly, VelocityGradientSlotsFollowHostOrder )
{
    GradientFixture fixture;
    TurbSaStateView state;

    FillSaVelocityGradient( fixture.Fields(), state );

    // component 0 = u, so its three derivatives occupy slots 0..2 in the
    // direction order x, y, z.
    EXPECT_EQ( state.velocityGradient[ 0 ], & fixture.storage[ 0 ][ 0 ] );
    EXPECT_EQ( state.velocityGradient[ 1 ], & fixture.storage[ 0 ][ 1 ] );
    EXPECT_EQ( state.velocityGradient[ 2 ], & fixture.storage[ 0 ][ 2 ] );
    EXPECT_EQ( state.velocityGradient[ 3 ], & fixture.storage[ 1 ][ 0 ] );
    EXPECT_EQ( state.velocityGradient[ 4 ], & fixture.storage[ 1 ][ 1 ] );
    EXPECT_EQ( state.velocityGradient[ 5 ], & fixture.storage[ 1 ][ 2 ] );
    EXPECT_EQ( state.velocityGradient[ 6 ], & fixture.storage[ 2 ][ 0 ] );
    EXPECT_EQ( state.velocityGradient[ 7 ], & fixture.storage[ 2 ][ 1 ] );
    EXPECT_EQ( state.velocityGradient[ 8 ], & fixture.storage[ 2 ][ 2 ] );
}

// The transposed mapping is the specific mistake to rule out: it would put
// dudx,dvdx,dwdx in slots 0..2.  Assert the slots are *not* that.
TEST( SaStateAssembly, VelocityGradientIsNotTransposed )
{
    GradientFixture fixture;
    TurbSaStateView state;

    FillSaVelocityGradient( fixture.Fields(), state );

    // A transpose would place the three velocities' x-derivatives first.
    const bool looksTransposed =
        state.velocityGradient[ 0 ] == & fixture.storage[ 0 ][ 0 ]
        && state.velocityGradient[ 1 ] == & fixture.storage[ 1 ][ 0 ]
        && state.velocityGradient[ 2 ] == & fixture.storage[ 2 ][ 0 ];
    EXPECT_FALSE( looksTransposed )
        << "slots 0..2 hold dudx,dvdx,dwdx: the mapping is transposed";

    // And state the positive form too, so the failure message is actionable.
    EXPECT_EQ( state.velocityGradient[ 1 ], & fixture.storage[ 0 ][ 1 ] )
        << "slot 1 must be dudy, i.e. gradient[u][y]";
    EXPECT_EQ( state.velocityGradient[ 3 ], & fixture.storage[ 1 ][ 0 ] )
        << "slot 3 must be dvdx, i.e. gradient[v][x]";
}

// Every slot must be filled.  A partially written view would read whatever the
// previous call left in the remaining slots.
TEST( SaStateAssembly, VelocityGradientFillsAllNineSlots )
{
    GradientFixture fixture;
    TurbSaStateView state;

    FillSaVelocityGradient( fixture.Fields(), state );

    for ( int slot = 0; slot < 9; ++ slot )
    {
        EXPECT_NE( state.velocityGradient[ slot ], nullptr )
            << "slot " << slot << " was not filled";
    }
}

// The sentinel values must be distinct, otherwise the transposition test above
// could not distinguish the two orders.  This guards the fixture itself.
TEST( SaStateAssembly, FixtureSentinelsAreDistinct )
{
    GradientFixture fixture;
    for ( int a = 0; a < 9; ++ a )
    {
        for ( int b = a + 1; b < 9; ++ b )
        {
            EXPECT_NE(
                fixture.storage[ a / 3 ][ a % 3 ],
                fixture.storage[ b / 3 ][ b % 3 ] )
                << "sentinels " << a << " and " << b << " collide";
        }
    }
}

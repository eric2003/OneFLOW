#include "ViscousFluxSelfCheck.h"

#include <gtest/gtest.h>

#include <limits>
#include <vector>

namespace
{

using namespace ONEFLOW;

// Five equation planes, three faces.  Values are chosen so that every face has
// a non-zero scale.
struct Planes
{
    std::vector< Real > storage[ 5 ];
    const Real * pointer[ 5 ];

    Planes()
    {
        for ( int equation = 0; equation < 5; ++ equation )
        {
            storage[ equation ].assign( 3, 0.0 );
        }
    }

    void Refresh()
    {
        for ( int equation = 0; equation < 5; ++ equation )
        {
            pointer[ equation ] = storage[ equation ].data();
        }
    }
};

TEST( ViscousFluxSelfCheck, IdenticalFluxPasses )
{
    Planes device;
    Planes host;
    for ( int equation = 0; equation < 5; ++ equation )
    {
        for ( int face = 0; face < 3; ++ face )
        {
            const Real value = 0.5 + equation + 2.0 * face;
            device.storage[ equation ][ face ] = value;
            host.storage[ equation ][ face ] = value;
        }
    }
    device.Refresh();
    host.Refresh();

    const ViscousFluxSelfCheckSummary summary = CompareViscousFaceFlux(
        device.pointer, host.pointer, 5, 3, 1e-10 );

    EXPECT_TRUE( summary.inputValid );
    EXPECT_TRUE( summary.Finite() );
    EXPECT_TRUE( summary.Passed() );
    EXPECT_EQ( summary.differingFaces, 0 );
    EXPECT_EQ( summary.comparedFaces, 3 );
    EXPECT_TRUE( summary.differences.empty() );
    EXPECT_EQ( summary.firstFace, -1 );
}

// Regression: the previous in-place check compared rel = |a-b|/scale against
// the tolerance.  Every comparison with NaN is false, so a NaN face looked like
// agreement.  Non-finite values must fail regardless of the tolerance.
TEST( ViscousFluxSelfCheck, NonFiniteDeviceValueFails )
{
    Planes device;
    Planes host;
    for ( int equation = 0; equation < 5; ++ equation )
    {
        for ( int face = 0; face < 3; ++ face )
        {
            device.storage[ equation ][ face ] = 1.0;
            host.storage[ equation ][ face ] = 1.0;
        }
    }
    device.storage[ 2 ][ 1 ] = std::numeric_limits< Real >::quiet_NaN();
    device.Refresh();
    host.Refresh();

    const ViscousFluxSelfCheckSummary summary = CompareViscousFaceFlux(
        device.pointer, host.pointer, 5, 3, 1e-10 );

    EXPECT_TRUE( summary.inputValid );
    EXPECT_FALSE( summary.Finite() );
    EXPECT_FALSE( summary.Passed() );
    EXPECT_EQ( summary.nonFiniteDeviceFaces, 1 );
    EXPECT_EQ( summary.nonFiniteHostFaces, 0 );
    EXPECT_EQ( summary.differingFaces, 1 );
    ASSERT_EQ( summary.differences.size(), 1u );
    EXPECT_EQ( summary.differences[ 0 ].face, 1 );
    EXPECT_TRUE( summary.differences[ 0 ].nonFinite );
}

TEST( ViscousFluxSelfCheck, NonFiniteHostValueFails )
{
    Planes device;
    Planes host;
    for ( int equation = 0; equation < 5; ++ equation )
    {
        for ( int face = 0; face < 3; ++ face )
        {
            device.storage[ equation ][ face ] = -2.0;
            host.storage[ equation ][ face ] = -2.0;
        }
    }
    host.storage[ 4 ][ 2 ] = std::numeric_limits< Real >::infinity();
    device.Refresh();
    host.Refresh();

    const ViscousFluxSelfCheckSummary summary = CompareViscousFaceFlux(
        device.pointer, host.pointer, 5, 3, 1e-10 );

    EXPECT_FALSE( summary.Passed() );
    EXPECT_EQ( summary.nonFiniteHostFaces, 1 );
    EXPECT_EQ( summary.nonFiniteDeviceFaces, 0 );
    EXPECT_EQ( summary.firstFace, 2 );
    ASSERT_EQ( summary.differences.size(), 1u );
    EXPECT_EQ( summary.differences[ 0 ].face, 2 );
}

// The threshold is relative to the largest magnitude on the face, not absolute.
TEST( ViscousFluxSelfCheck, ToleranceIsRelativeToFaceScale )
{
    Planes device;
    Planes host;
    for ( int equation = 0; equation < 5; ++ equation )
    {
        device.storage[ equation ][ 0 ] = 1.0 * ( equation + 1 );
        host.storage[ equation ][ 0 ] = device.storage[ equation ][ 0 ];
        device.storage[ equation ][ 1 ] = 1000.0 * ( equation + 1 );
        host.storage[ equation ][ 1 ] = device.storage[ equation ][ 1 ];
        device.storage[ equation ][ 2 ] = 1.0 * ( equation + 1 );
        host.storage[ equation ][ 2 ] = device.storage[ equation ][ 2 ];
    }
    // Face 0 has scale 5: 1e-10 is rel = 2e-11, accepted.
    device.storage[ 4 ][ 0 ] += 1e-10;
    // The same absolute difference is applied to the small face (scale 5) and
    // the large face (scale 5000).  Only the small face exceeds 1e-10.
    device.storage[ 4 ][ 1 ] += 1.2e-9; // rel = 2.4e-13
    device.storage[ 4 ][ 2 ] += 1.2e-9; // rel = 2.4e-10
    device.Refresh();
    host.Refresh();

    const ViscousFluxSelfCheckSummary summary = CompareViscousFaceFlux(
        device.pointer, host.pointer, 5, 3, 1e-10 );

    EXPECT_FALSE( summary.Passed() );
    EXPECT_EQ( summary.differingFaces, 1 );
    EXPECT_EQ( summary.firstFace, 2 );
    EXPECT_EQ( summary.firstEquation, 4 );
    EXPECT_EQ( summary.worstFace, 2 );
    EXPECT_EQ( summary.worstEquation, 4 );
    EXPECT_NEAR( summary.worstRelative, 2.4e-10, 1e-15 );
}

TEST( ViscousFluxSelfCheck, AllZeroFacesAreNotCompared )
{
    Planes device;
    Planes host;
    device.Refresh();
    host.Refresh();

    const ViscousFluxSelfCheckSummary summary = CompareViscousFaceFlux(
        device.pointer, host.pointer, 5, 3, 1e-10 );

    EXPECT_TRUE( summary.Passed() );
    EXPECT_EQ( summary.comparedFaces, 0 );
    EXPECT_EQ( summary.differingFaces, 0 );
}

// A face whose device value is exactly zero and whose host value is tiny still
// differs: the scale is the larger magnitude, so relative error is 1.
TEST( ViscousFluxSelfCheck, ZeroVersusTinyValueIsADifference )
{
    Planes device;
    Planes host;
    device.Refresh();
    host.Refresh();
    host.storage[ 0 ][ 0 ] = 1e-30;
    host.Refresh();

    const ViscousFluxSelfCheckSummary summary = CompareViscousFaceFlux(
        device.pointer, host.pointer, 5, 3, 1e-10 );

    EXPECT_FALSE( summary.Passed() );
    EXPECT_EQ( summary.differingFaces, 1 );
    EXPECT_DOUBLE_EQ( summary.worstRelative, 1.0 );
}

TEST( ViscousFluxSelfCheck, ReportsFirstAndWorstDifference )
{
    Planes device;
    Planes host;
    for ( int equation = 0; equation < 5; ++ equation )
    {
        for ( int face = 0; face < 3; ++ face )
        {
            device.storage[ equation ][ face ] = 1.0;
            host.storage[ equation ][ face ] = 1.0;
        }
    }
    // Face 0 differs slightly in equation 1.
    device.storage[ 1 ][ 0 ] = 1.0 + 1e-6;
    // Face 2 differs most, in equation 3.
    device.storage[ 3 ][ 2 ] = 1.0 - 1e-2;
    device.Refresh();
    host.Refresh();

    const ViscousFluxSelfCheckSummary summary = CompareViscousFaceFlux(
        device.pointer, host.pointer, 5, 3, 1e-10 );

    EXPECT_FALSE( summary.Passed() );
    EXPECT_EQ( summary.differingFaces, 2 );
    EXPECT_EQ( summary.firstFace, 0 );
    EXPECT_EQ( summary.firstEquation, 1 );
    EXPECT_EQ( summary.worstFace, 2 );
    EXPECT_EQ( summary.worstEquation, 3 );
    EXPECT_NEAR( summary.worstRelative, 1e-2, 1e-13 );
    ASSERT_EQ( summary.differences.size(), 2u );
    EXPECT_EQ( summary.differences[ 0 ].face, 0 );
    EXPECT_EQ( summary.differences[ 1 ].face, 2 );
}

TEST( ViscousFluxSelfCheck, DifferenceRecordingRespectsTheRequestedCap )
{
    Planes device;
    Planes host;
    for ( int equation = 0; equation < 5; ++ equation )
    {
        device.storage[ equation ].assign( 10, 1.0 );
        host.storage[ equation ].assign( 10, 1.0 );
    }
    for ( int face = 0; face < 10; ++ face )
    {
        device.storage[ 0 ][ face ] = 2.0;
    }
    device.Refresh();
    host.Refresh();

    const ViscousFluxSelfCheckSummary summary = CompareViscousFaceFlux(
        device.pointer, host.pointer, 5, 10, 1e-10, 2 );

    EXPECT_EQ( summary.differingFaces, 10 );
    ASSERT_EQ( summary.differences.size(), 2u );
    EXPECT_EQ( summary.differences[ 0 ].face, 0 );
    EXPECT_EQ( summary.differences[ 1 ].face, 1 );
}

TEST( ViscousFluxSelfCheck, RejectsUnusableInput )
{
    Planes device;
    Planes host;
    device.Refresh();
    host.Refresh();

    const Real nan = std::numeric_limits< Real >::quiet_NaN();

    EXPECT_FALSE( CompareViscousFaceFlux(
        nullptr, host.pointer, 5, 3, 1e-10 ).Passed() );
    EXPECT_FALSE( CompareViscousFaceFlux(
        device.pointer, nullptr, 5, 3, 1e-10 ).Passed() );

    const Real * missingPlane[ 5 ] = {
        device.pointer[ 0 ], device.pointer[ 1 ], nullptr,
        device.pointer[ 3 ], device.pointer[ 4 ] };
    const ViscousFluxSelfCheckSummary missing = CompareViscousFaceFlux(
        missingPlane, host.pointer, 5, 3, 1e-10 );
    EXPECT_FALSE( missing.inputValid );
    EXPECT_FALSE( missing.Passed() );

    EXPECT_FALSE( CompareViscousFaceFlux(
        device.pointer, host.pointer, 0, 3, 1e-10 ).Passed() );
    EXPECT_FALSE( CompareViscousFaceFlux(
        device.pointer, host.pointer, 5, -1, 1e-10 ).Passed() );
    EXPECT_FALSE( CompareViscousFaceFlux(
        device.pointer, host.pointer, 5, 3, nan ).Passed() );
    EXPECT_FALSE( CompareViscousFaceFlux(
        device.pointer, host.pointer, 5, 3, -1e-10 ).Passed() );
}

TEST( ViscousFluxSelfCheck, RecognisesNonFiniteValues )
{
    EXPECT_TRUE( IsViscousFluxValueFinite( 0.0 ) );
    EXPECT_TRUE( IsViscousFluxValueFinite( -3.5 ) );
    EXPECT_FALSE( IsViscousFluxValueFinite(
        std::numeric_limits< Real >::quiet_NaN() ) );
    EXPECT_FALSE( IsViscousFluxValueFinite(
        std::numeric_limits< Real >::infinity() ) );
    EXPECT_FALSE( IsViscousFluxValueFinite(
        -std::numeric_limits< Real >::infinity() ) );
}

}

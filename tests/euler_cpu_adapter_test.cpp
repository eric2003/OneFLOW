#include "EulerCpuAdapter.h"

#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <vector>

namespace
{

using namespace ONEFLOW;

TEST( EulerCpuAdapter, ConvertsThreeEquationPrimitiveState )
{
    const Real primitive[] = { 1.2, 0.25, 0.9 };
    Real conserved[ 3 ] = {};
    PrimitiveToConserved( primitive, 3, 1.4, conserved );

    EXPECT_DOUBLE_EQ( conserved[ 0 ], 1.2 );
    EXPECT_DOUBLE_EQ( conserved[ 1 ], 0.3 );
    EXPECT_DOUBLE_EQ( conserved[ 2 ], 0.9 / 0.4 + 0.5 * 1.2 * 0.25 * 0.25 );
}

TEST( EulerCpuAdapter, PacksFiveEquationFaceBatchEquationMajor )
{
    constexpr int nFaces = 2;
    const Real primitiveLeft[] = {
        1.0, 0.8, 0.2, 0.4, 0.3, 0.1, 0.0, 0.2, 1.0, 0.9 };
    const Real primitiveRight[] = {
        0.9, 1.1, 0.1, 0.5, 0.0, 0.2, 0.2, 0.0, 0.8, 1.2 };
    Real conservedLeft[ 5 * nFaces ] = {};
    Real conservedRight[ 5 * nFaces ] = {};

    PackPrimitiveFaceStates( primitiveLeft, primitiveRight, nFaces, 5, 1.4,
        conservedLeft, conservedRight );

    EXPECT_DOUBLE_EQ( conservedLeft[ 0 ], 1.0 );
    EXPECT_DOUBLE_EQ( conservedLeft[ 1 ], 0.8 );
    EXPECT_DOUBLE_EQ( conservedLeft[ 2 ], 0.2 );
    EXPECT_DOUBLE_EQ( conservedLeft[ 3 ], 0.32 );
    EXPECT_DOUBLE_EQ( conservedLeft[ 4 ], 0.3 );
    EXPECT_DOUBLE_EQ( conservedLeft[ 5 ], 0.08 );
    EXPECT_DOUBLE_EQ( conservedLeft[ 6 ], 0.0 );
    EXPECT_DOUBLE_EQ( conservedLeft[ 7 ], 0.16 );
    EXPECT_GT( conservedLeft[ 8 ], 0.0 );
    EXPECT_GT( conservedLeft[ 9 ], 0.0 );
    EXPECT_NE( conservedLeft[ 8 ], conservedRight[ 8 ] );
}

TEST( EulerCpuAdapter, RejectsNonPhysicalOrUnsupportedState )
{
    const Real negativeDensity[] = { -1.0, 0.0, 1.0 };
    Real conserved[ 5 ] = {};
    EXPECT_THROW( PrimitiveToConserved( negativeDensity, 3, 1.4, conserved ),
        std::invalid_argument );

    const Real unsupported[] = { 1.0, 0.0, 0.0, 1.0 };
    EXPECT_THROW( PrimitiveToConserved( unsupported, 4, 1.4, conserved ),
        std::invalid_argument );
}

TEST( EulerCpuAdapter, MapsResidualUsingExplicitBoundaryMask )
{
    constexpr int nFaces = 3;
    constexpr int nCells = 5;
    Real faceValues[] = { 1.0, 2.0, 3.0, 10.0, 20.0, 30.0, 100.0, 200.0, 300.0 };
    int left[] = { 0, 1, 2 };
    int right[] = { 2, 3, 4 };
    unsigned char boundaryMask[] = { 0, 1, 0 };
    Real residualValues[ 3 * nCells ] = {};

    FaceFluxView flux{ nFaces, 3, faceValues };
    FaceConnectivityView connectivity;
    connectivity.nFaces = nFaces;
    connectivity.nBoundaryFaces = 1;
    connectivity.leftCell = left;
    connectivity.rightCell = right;
    connectivity.boundaryMask = boundaryMask;
    ResidualView residual{ nCells, 3, residualValues };

    EulerCpuAdapter adapter;
    EXPECT_NO_THROW( adapter.AddFaceFlux( flux, connectivity, residual ) );
    EXPECT_DOUBLE_EQ( residualValues[ 0 ], -1.0 );
    EXPECT_DOUBLE_EQ( residualValues[ 1 ], -2.0 );
    EXPECT_DOUBLE_EQ( residualValues[ 2 ], -2.0 );
    EXPECT_DOUBLE_EQ( residualValues[ 4 ], 3.0 );
    EXPECT_DOUBLE_EQ( residualValues[ 5 ], -10.0 );
    EXPECT_DOUBLE_EQ( residualValues[ 7 ], -20.0 );
    EXPECT_DOUBLE_EQ( residualValues[ 9 ], 30.0 );
}

TEST( EulerCpuAdapter, ComputesBatchFluxAfterConversion )
{
    constexpr int nFaces = 2;
    const Real primitiveLeft[] = { 1.0, 1.0, 0.25, 0.25, 1.0, 1.0 };
    const Real primitiveRight[] = { 1.0, 1.0, 0.25, 0.25, 1.0, 1.0 };
    const Real normal[] = { 1.0, 1.0 };
    const Real area[] = { 2.0, 3.0 };
    Real values[ 3 * nFaces ] = {};

    PrimitiveFaceStateView state;
    state.nFaces = nFaces;
    state.nEquations = 3;
    state.primitiveLeft = primitiveLeft;
    state.primitiveRight = primitiveRight;
    state.xNormal = normal;
    state.faceArea = area;

    FaceFluxView flux{ nFaces, 3, values };
    EulerCpuAdapter adapter;
    EXPECT_NO_THROW( adapter.CalcInvFlux( state, flux ) );
    for ( Real value : values ) EXPECT_TRUE( std::isfinite( value ) );
    EXPECT_GT( values[ 0 ], 0.0 );
    EXPECT_GT( values[ 1 ], values[ 0 ] );
}


TEST( EulerCpuAdapter, ChecksPhysicalityAndInternalFaceConservation )
{
    constexpr int nFaces = 2;
    constexpr int nCells = 3;
    constexpr int nEquations = 5;
    const Real primitiveLeft[] = {
        1.0, 1.1, 0.2, -0.1, 1.0,
        0.9, 0.8, -0.3, 0.4, 0.7 };
    const Real primitiveRight[] = {
        0.95, 1.0, 0.1, -0.2, 0.9,
        1.05, 0.7, -0.1, 0.3, 0.8 };
    const Real normalX[] = { 1.0, 1.0 };
    const Real normalY[] = { 0.0, 0.0 };
    const Real normalZ[] = { 0.0, 0.0 };
    const Real area[] = { 1.0, 2.0 };
    Real fluxValues[ nFaces * nEquations ] = {};

    PrimitiveFaceStateView state;
    state.nFaces = nFaces;
    state.nEquations = nEquations;
    state.primitiveLeft = primitiveLeft;
    state.primitiveRight = primitiveRight;
    state.xNormal = normalX;
    state.yNormal = normalY;
    state.zNormal = normalZ;
    state.faceArea = area;
    state.gamma = 1.4;

    FaceFluxView flux{ nFaces, nEquations, fluxValues };
    EulerCpuAdapter adapter;
    ASSERT_NO_THROW( adapter.CalcInvFlux( state, flux, 1 ) );
    for ( Real value : fluxValues )
    {
        EXPECT_TRUE( std::isfinite( value ) );
    }

    int left[] = { 0, 1 };
    int right[] = { 1, 2 };
    unsigned char boundaryMask[] = { 0, 0 };
    Real residualValues[ nCells * nEquations ] = {};
    FaceConnectivityView connectivity{
        nFaces, 0, left, right, boundaryMask };
    ResidualView residual{ nCells, nEquations, residualValues };
    ASSERT_NO_THROW( adapter.AddFaceFlux( flux, connectivity, residual ) );

    for ( int equation = 0; equation < nEquations; ++ equation )
    {
        Real sum = 0.0;
        for ( int cell = 0; cell < nCells; ++ cell )
        {
            sum += residualValues[ equation * nCells + cell ];
        }
        EXPECT_NEAR( sum, 0.0, 1.0e-14 );
    }
}

TEST( EulerCpuAdapter, MatchesOneflowLaxFriedrichsReference )
{
    const Real primitiveLeft[] = { 1.0, 2.0, 0.3, -0.2, 1.0 };
    const Real primitiveRight[] = { 0.8, 1.0, -0.1, 0.4, 0.7 };
    const Real normal[] = { 0.6 };
    const Real tangent[] = { 0.8 };
    const Real zero[] = { 0.0 };
    const Real meshVelocity[] = { 0.05 };
    const Real area[] = { 2.0 };
    Real values[ 5 ] = {};

    PrimitiveFaceStateView state;
    state.nFaces = 1;
    state.nEquations = 5;
    state.primitiveLeft = primitiveLeft;
    state.primitiveRight = primitiveRight;
    state.xNormal = normal;
    state.yNormal = tangent;
    state.zNormal = zero;
    state.meshVelocityNormal = meshVelocity;
    state.faceArea = area;
    state.gamma = 1.4;

    FaceFluxView flux{ 1, 5, values };
    EulerCpuAdapter adapter;
    ASSERT_NO_THROW( adapter.CalcInvFlux( state, flux, 1 ) );

    const Real expected[] = {
        2.1889497342723656,
        6.713698405634195,
        2.5430044951174953,
        -1.227269309108151,
        14.155125131686214
    };
    for ( int equation = 0; equation < 5; ++ equation )
    {
        EXPECT_NEAR( values[ equation ], expected[ equation ], 1.0e-12 );
    }
}

} // namespace

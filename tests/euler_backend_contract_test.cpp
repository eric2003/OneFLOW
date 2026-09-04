#include "OneDEulerBackend.h"

#include <gtest/gtest.h>

#ifdef ONEFLOW_1D_USE_HIP
#include <hip/hip_runtime.h>
#endif

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

using oneflow_1d::CpuEulerBackend;
using oneflow_1d::EulerBackend;
using oneflow_1d::EulerBoundary;
using oneflow_1d::EulerComponents;
using oneflow_1d::EulerProblem;
using oneflow_1d::EulerRunMode;
using oneflow_1d::EulerRunOptions;
using oneflow_1d::EulerTrace;

constexpr double kGamma = 1.4;
constexpr int kNx = 32;
constexpr double kDt = 0.0005;
constexpr double kDx = 1.0 / kNx;
constexpr double kPi = 3.14159265358979323846;

EulerProblem Problem()
{
    return { kNx, kGamma, kDt, kDx, EulerBoundary::Periodic };
}

std::vector< double > InitialState()
{
    std::vector< double > state( EulerComponents * kNx );
    for ( int cell = 0; cell < kNx; ++ cell )
    {
        const double x = ( cell + 0.5 ) * kDx;
        const double density = 1.0 + 0.08 * std::sin( 2.0 * kPi * x );
        const double velocity = 0.2 + 0.04 * std::cos( 2.0 * kPi * x );
        const double pressure = 1.0 + 0.05 * std::sin( 4.0 * kPi * x );
        state[ cell ] = density;
        state[ kNx + cell ] = density * velocity;
        state[ 2 * kNx + cell ] = pressure / ( kGamma - 1.0 )
            + 0.5 * density * velocity * velocity;
    }
    return state;
}

std::unique_ptr< EulerBackend > MakeBackend()
{
#ifdef ONEFLOW_1D_USE_HIP
    return std::make_unique< oneflow_1d::HipEulerBackend >();
#else
    return std::make_unique< CpuEulerBackend >();
#endif
}

std::vector< double > RunNoTrace(
    const EulerBackend & backend, const std::vector< double > & initial,
    int steps )
{
    auto state = backend.CreateState( Problem() );
    backend.Upload( *state, initial.data() );
    backend.Advance( *state, steps, {} );
    std::vector< double > result( initial.size() );
    backend.Download( *state, result.data() );
    return result;
}

void ExpectClose(
    const std::vector< double > & expected,
    const std::vector< double > & actual, const std::string & label )
{
    ASSERT_EQ( expected.size(), actual.size() );
    for ( std::size_t index = 0; index < expected.size(); ++ index )
    {
        const double reference = expected[ index ];
        const double tolerance = 1.0e-15
            + 1.0e-15 * std::max( 1.0, std::abs( reference ) );
        EXPECT_LE( std::abs( actual[ index ] - reference ), tolerance )
            << label << " at index " << index << ", expected " << reference
            << ", actual " << actual[ index ];
    }
}

void ExpectPhysical( const std::vector< double > & state )
{
    ASSERT_EQ( state.size(), static_cast< std::size_t >( 3 * kNx ) );
    for ( int cell = 0; cell < kNx; ++ cell )
    {
        const double density = state[ cell ];
        const double momentum = state[ kNx + cell ];
        const double energy = state[ 2 * kNx + cell ];
        const double pressure = ( kGamma - 1.0 )
            * ( energy - 0.5 * momentum * momentum / density );
        EXPECT_TRUE( std::isfinite( density ) );
        EXPECT_TRUE( std::isfinite( momentum ) );
        EXPECT_TRUE( std::isfinite( energy ) );
        EXPECT_GT( density, 0.0 );
        EXPECT_GT( pressure, 0.0 );
    }
}

TEST( EulerBackendContract, ReportsBackendIdentity )
{
    const auto backend = MakeBackend();
#ifdef ONEFLOW_1D_USE_HIP
    EXPECT_STREQ( backend->Name(), "HIP" );
    EXPECT_TRUE( backend->IsAccelerator() );
#else
    EXPECT_STREQ( backend->Name(), "CPU" );
    EXPECT_FALSE( backend->IsAccelerator() );
#endif
}

#ifdef ONEFLOW_1D_USE_HIP
TEST( EulerBackendContract, RequiresVisibleHipDevice )
{
    int deviceCount = 0;
    ASSERT_EQ( hipGetDeviceCount( &deviceCount ), hipSuccess );
    ASSERT_GT( deviceCount, 0 );
}
#endif

TEST( EulerBackendContract, FullTraceMatchesCpuReference )
{
    const auto backend = MakeBackend();
    const std::vector< double > initial = InitialState();
    EulerTrace actual;
    backend->Step(
        initial.data(), kNx, kGamma, kDt, kDx, EulerBoundary::Periodic, actual );

    CpuEulerBackend referenceBackend;
    EulerTrace expected;
    referenceBackend.Step(
        initial.data(), kNx, kGamma, kDt, kDx, EulerBoundary::Periodic, expected );

    ASSERT_EQ( actual.nx, kNx );
    ASSERT_EQ( actual.faceLeft.size(), expected.faceLeft.size() );
    ASSERT_EQ( actual.faceRight.size(), expected.faceRight.size() );
    ASSERT_EQ( actual.numericalFlux.size(), expected.numericalFlux.size() );
    ASSERT_EQ( actual.residual.size(), expected.residual.size() );
    ASSERT_EQ( actual.state.size(), expected.state.size() );
    ExpectClose( expected.faceLeft, actual.faceLeft, "faceLeft" );
    ExpectClose( expected.faceRight, actual.faceRight, "faceRight" );
    ExpectClose( expected.numericalFlux, actual.numericalFlux, "flux" );
    ExpectClose( expected.residual, actual.residual, "residual" );
    ExpectClose( expected.state, actual.state, "state" );
    ExpectPhysical( std::vector< double >(
        actual.state.end() - 3 * kNx, actual.state.end() ) );
}

TEST( EulerBackendContract, NoTraceMatchesCpuReference )
{
    const auto backend = MakeBackend();
    const std::vector< double > initial = InitialState();
    const std::vector< double > actual = RunNoTrace( *backend, initial, 3 );
    const std::vector< double > expected = RunNoTrace(
        CpuEulerBackend(), initial, 3 );
    ExpectClose( expected, actual, "final state" );
    ExpectPhysical( actual );
}

TEST( EulerBackendContract, LifecycleCanBeReused )
{
    const auto backend = MakeBackend();
    const std::vector< double > initial = InitialState();
    auto state = backend->CreateState( Problem() );
    backend->Upload( *state, initial.data() );

    backend->Advance( *state, 1, {} );
    std::vector< double > first( initial.size() );
    backend->Download( *state, first.data() );
    backend->Advance( *state, 1, {} );
    std::vector< double > second( initial.size() );
    backend->Download( *state, second.data() );

    ExpectClose( RunNoTrace( CpuEulerBackend(), initial, 1 ), first, "first" );
    ExpectClose( RunNoTrace( CpuEulerBackend(), initial, 2 ), second, "second" );
}

TEST( EulerBackendContract, RejectsInvalidRequests )
{
    const auto backend = MakeBackend();
    EulerProblem invalid = Problem();
    invalid.nx = 0;
    EXPECT_THROW( backend->CreateState( invalid ), std::invalid_argument );

    auto state = backend->CreateState( Problem() );
    EXPECT_THROW( backend->Upload( *state, nullptr ), std::invalid_argument );
    EXPECT_THROW( backend->Advance( *state, 1, {} ), std::logic_error );

    const std::vector< double > initial = InitialState();
    backend->Upload( *state, initial.data() );
    EXPECT_THROW( backend->Advance( *state, 0, {} ), std::invalid_argument );

    EulerTrace trace;
    EulerRunOptions fullTrace = { EulerRunMode::FullTrace, &trace, nullptr };
    EXPECT_THROW(
        backend->Advance( *state, 2, fullTrace ), std::invalid_argument );
    EulerRunOptions noTrace = { EulerRunMode::NoTrace, &trace, nullptr };
    EXPECT_THROW(
        backend->Advance( *state, 1, noTrace ), std::invalid_argument );
}

} // namespace

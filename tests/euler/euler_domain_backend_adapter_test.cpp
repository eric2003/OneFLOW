#include "OneDEulerDomainAdapter.h"

#include <gtest/gtest.h>

#ifdef ONEFLOW_1D_USE_HIP
#include "AccelRuntime.h"
#endif

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace
{

constexpr int kCells = 32;
constexpr double kGamma = 1.4;
constexpr double kDt = 0.0005;
constexpr double kDx = 1.0 / kCells;

ONEFLOW::EulerDomainProblem SharedProblem()
{
    ONEFLOW::EulerDomainProblem problem;
    problem.nCells = kCells;
    problem.nGhostCells = 0;
    problem.nEquations = oneflow_1d::EulerComponents;
    problem.gamma = kGamma;
    problem.dt = kDt;
    problem.dx = kDx;
    problem.boundary = ONEFLOW::EulerDomainBoundary::Periodic;
    return problem;
}

ONEFLOW::EulerDomainStateKey CpuKey()
{
    return { 0, 0, 0, ONEFLOW::AccelBackendKind::CPU, -1 };
}

#ifdef ONEFLOW_1D_USE_HIP
ONEFLOW::EulerDomainStateKey HipKey( int deviceId = 0 )
{
    return { 0, 0, 0, ONEFLOW::AccelBackendKind::HIP, deviceId };
}

void EnsureHipRuntime()
{
    static const bool initialized = []() {
        ONEFLOW::InitializeAccelRuntime( 0, 1 );
        return true;
    }();
    static_cast< void >( initialized );
}
#endif

std::vector< double > InitialState()
{
    std::vector< double > values(
        oneflow_1d::EulerComponents * kCells );
    for ( int cell = 0; cell < kCells; ++ cell )
    {
        const double x = ( cell + 0.5 ) * kDx;
        const double density = 1.0 + 0.05 * std::sin( 6.283185307179586 * x );
        const double velocity = 0.2;
        const double pressure = 1.0;
        values[ cell ] = density;
        values[ kCells + cell ] = density * velocity;
        values[ 2 * kCells + cell ] =
            pressure / ( kGamma - 1.0 )
            + 0.5 * density * velocity * velocity;
    }
    return values;
}

void ExpectClose(
    const std::vector< double > & expected,
    const std::vector< double > & actual,
    const std::string & label )
{
    ASSERT_EQ( expected.size(), actual.size() );
    for ( std::size_t index = 0; index < expected.size(); ++ index )
    {
        const double reference = expected[ index ];
        const double tolerance =
            1.0e-15 + 1.0e-15 * std::max( 1.0, std::abs( reference ) );
        EXPECT_LE( std::abs( actual[ index ] - reference ), tolerance )
            << label << " at index " << index;
    }
}

void ExpectTraceClose(
    const oneflow_1d::EulerTrace & expected,
    const oneflow_1d::EulerTrace & actual )
{
    EXPECT_EQ( actual.nx, expected.nx );
    ExpectClose( expected.faceLeft, actual.faceLeft, "faceLeft" );
    ExpectClose( expected.faceRight, actual.faceRight, "faceRight" );
    ExpectClose( expected.numericalFlux, actual.numericalFlux, "flux" );
    ExpectClose( expected.residual, actual.residual, "residual" );
    ExpectClose( expected.state, actual.state, "state" );
}

TEST( EulerDomainBackendAdapter, MatchesDirectCpuLifecycle )
{
    oneflow_1d::CpuEulerBackend direct;
    oneflow_1d::EulerDomainBackendAdapter adapter(
        direct, ONEFLOW::AccelBackendKind::CPU, -1 );
    const auto sharedProblem = SharedProblem();
    const auto initial = InitialState();

    auto sharedState = adapter.CreateState( sharedProblem, CpuKey() );
    ONEFLOW::EulerDomainConstFieldView input{
        kCells, oneflow_1d::EulerComponents, initial.data() };
    adapter.Upload( *sharedState, input );
    adapter.Advance( *sharedState, 2, {} );
    std::vector< double > sharedFinal( initial.size() );
    ONEFLOW::EulerDomainFieldView output{
        kCells, oneflow_1d::EulerComponents, sharedFinal.data() };
    adapter.Download( *sharedState, output );

    const oneflow_1d::EulerProblem directProblem{
        kCells, kGamma, kDt, kDx,
        oneflow_1d::EulerBoundary::Periodic,
        oneflow_1d::EulerMethod::Rusanov };
    auto directState = direct.CreateState( directProblem );
    direct.Upload( *directState, initial.data() );
    direct.Advance( *directState, 2, {} );
    std::vector< double > directFinal( initial.size() );
    direct.Download( *directState, directFinal.data() );

    EXPECT_EQ( sharedFinal, directFinal );
    EXPECT_STREQ( adapter.Name(), direct.Name() );
    EXPECT_FALSE( adapter.IsAccelerator() );
}

TEST( EulerDomainBackendAdapter, RejectsUnsupportedProblemAndOwnership )
{
    oneflow_1d::CpuEulerBackend cpu;
    oneflow_1d::EulerDomainBackendAdapter adapter(
        cpu, ONEFLOW::AccelBackendKind::CPU, -1 );

    auto problem = SharedProblem();
    problem.nEquations = 5;
    EXPECT_THROW(
        adapter.CreateState( problem, CpuKey() ), std::invalid_argument );

    problem = SharedProblem();
    problem.nGhostCells = 1;
    EXPECT_THROW(
        adapter.CreateState( problem, CpuKey() ), std::invalid_argument );

    problem = SharedProblem();
    problem.boundary = ONEFLOW::EulerDomainBoundary::Wall;
    EXPECT_THROW(
        adapter.CreateState( problem, CpuKey() ), std::invalid_argument );

    problem = SharedProblem();
    auto key = CpuKey();
    key.backend = ONEFLOW::AccelBackendKind::HIP;
    key.deviceId = 0;
    EXPECT_THROW(
        adapter.CreateState( problem, key ), std::invalid_argument );
}

TEST( EulerDomainBackendAdapter, RejectsExternalStageScheduling )
{
    oneflow_1d::CpuEulerBackend cpu;
    oneflow_1d::EulerDomainBackendAdapter adapter(
        cpu, ONEFLOW::AccelBackendKind::CPU, -1 );
    const auto problem = SharedProblem();
    auto state = adapter.CreateState( problem, CpuKey() );
    const auto initial = InitialState();
    ONEFLOW::EulerDomainConstFieldView input{
        kCells, oneflow_1d::EulerComponents, initial.data() };
    adapter.Upload( *state, input );

    ONEFLOW::EulerDomainRunOptions options;
    options.stageCount = 3;
    EXPECT_THROW(
        adapter.Advance( *state, 1, options ), std::invalid_argument );
}

#ifdef ONEFLOW_1D_USE_HIP
TEST( EulerDomainBackendAdapter, HipNoTraceMatchesDirectLifecycleAndStats )
{
    EnsureHipRuntime();
    oneflow_1d::HipEulerBackend hip;
    oneflow_1d::EulerDomainBackendAdapter adapter(
        hip, ONEFLOW::AccelBackendKind::HIP, 0 );
    const auto problem = SharedProblem();
    const auto initial = InitialState();

    auto sharedState = adapter.CreateState( problem, HipKey() );
    ONEFLOW::EulerDomainConstFieldView input{
        kCells, oneflow_1d::EulerComponents, initial.data() };
    adapter.Upload( *sharedState, input );
    oneflow_1d::EulerRunStats sharedStats;
    ONEFLOW::EulerDomainRunOptions sharedOptions;
    sharedOptions.stats = &sharedStats;
    adapter.Advance( *sharedState, 2, sharedOptions );
    std::vector< double > sharedFinal( initial.size() );
    ONEFLOW::EulerDomainFieldView output{
        kCells, oneflow_1d::EulerComponents, sharedFinal.data() };
    adapter.Download( *sharedState, output );

    const oneflow_1d::EulerProblem directProblem{
        kCells, kGamma, kDt, kDx,
        oneflow_1d::EulerBoundary::Periodic,
        oneflow_1d::EulerMethod::Rusanov };
    auto directState = hip.CreateState( directProblem );
    hip.Upload( *directState, initial.data() );
    oneflow_1d::EulerRunStats directStats;
    hip.Advance(
        *directState, 2,
        { oneflow_1d::EulerRunMode::NoTrace, nullptr, &directStats } );
    std::vector< double > directFinal( initial.size() );
    hip.Download( *directState, directFinal.data() );

    ExpectClose( directFinal, sharedFinal, "HIP NoTrace final state" );
    EXPECT_EQ( sharedStats.kernelLaunches, directStats.kernelLaunches );
    EXPECT_EQ(
        sharedStats.synchronizationCount,
        directStats.synchronizationCount );
    EXPECT_EQ( sharedStats.kernelLaunches, 12 );
    EXPECT_EQ( sharedStats.synchronizationCount, 1 );
    EXPECT_GE( sharedStats.kernelMilliseconds, 0.0 );
    EXPECT_STREQ( adapter.Name(), hip.Name() );
    EXPECT_TRUE( adapter.IsAccelerator() );
}

TEST( EulerDomainBackendAdapter, HipFullTraceMatchesDirectPath )
{
    EnsureHipRuntime();
    oneflow_1d::HipEulerBackend hip;
    oneflow_1d::EulerDomainBackendAdapter adapter(
        hip, ONEFLOW::AccelBackendKind::HIP, 0 );
    const auto problem = SharedProblem();
    const auto initial = InitialState();

    auto sharedState = adapter.CreateState( problem, HipKey() );
    ONEFLOW::EulerDomainConstFieldView input{
        kCells, oneflow_1d::EulerComponents, initial.data() };
    adapter.Upload( *sharedState, input );
    oneflow_1d::EulerTrace sharedTrace;
    ONEFLOW::EulerDomainRunOptions sharedOptions;
    sharedOptions.mode = ONEFLOW::EulerDomainRunMode::FullTrace;
    sharedOptions.trace = &sharedTrace;
    adapter.Advance( *sharedState, 1, sharedOptions );

    const oneflow_1d::EulerProblem directProblem{
        kCells, kGamma, kDt, kDx,
        oneflow_1d::EulerBoundary::Periodic,
        oneflow_1d::EulerMethod::Rusanov };
    auto directState = hip.CreateState( directProblem );
    hip.Upload( *directState, initial.data() );
    oneflow_1d::EulerTrace directTrace;
    hip.Advance(
        *directState, 1,
        { oneflow_1d::EulerRunMode::FullTrace, &directTrace, nullptr } );

    ExpectTraceClose( directTrace, sharedTrace );
}

TEST( EulerDomainBackendAdapter, HipEnforcesDeviceIdentityAndLifecycleReuse )
{
    EnsureHipRuntime();
    oneflow_1d::HipEulerBackend hip;
    oneflow_1d::EulerDomainBackendAdapter adapter(
        hip, ONEFLOW::AccelBackendKind::HIP, 0 );
    oneflow_1d::EulerDomainBackendAdapter otherDevice(
        hip, ONEFLOW::AccelBackendKind::HIP, 1 );
    const auto problem = SharedProblem();
    const auto initial = InitialState();

    EXPECT_THROW(
        adapter.CreateState( problem, HipKey( 1 ) ),
        std::invalid_argument );

    auto state = adapter.CreateState( problem, HipKey() );
    ONEFLOW::EulerDomainConstFieldView input{
        kCells, oneflow_1d::EulerComponents, initial.data() };
    adapter.Upload( *state, input );
    EXPECT_THROW(
        otherDevice.Advance( *state, 1, {} ),
        std::invalid_argument );

    adapter.Advance( *state, 1, {} );
    std::vector< double > first( initial.size() );
    ONEFLOW::EulerDomainFieldView firstOutput{
        kCells, oneflow_1d::EulerComponents, first.data() };
    adapter.Download( *state, firstOutput );
    adapter.Advance( *state, 1, {} );
    std::vector< double > second( initial.size() );
    ONEFLOW::EulerDomainFieldView secondOutput{
        kCells, oneflow_1d::EulerComponents, second.data() };
    adapter.Download( *state, secondOutput );

    const oneflow_1d::EulerProblem directProblem{
        kCells, kGamma, kDt, kDx,
        oneflow_1d::EulerBoundary::Periodic,
        oneflow_1d::EulerMethod::Rusanov };
    auto directState = hip.CreateState( directProblem );
    hip.Upload( *directState, initial.data() );
    hip.Advance( *directState, 1, {} );
    std::vector< double > directFirst( initial.size() );
    hip.Download( *directState, directFirst.data() );
    hip.Advance( *directState, 1, {} );
    std::vector< double > directSecond( initial.size() );
    hip.Download( *directState, directSecond.data() );

    ExpectClose( directFirst, first, "HIP first lifecycle advance" );
    ExpectClose( directSecond, second, "HIP second lifecycle advance" );
}
#endif

} // namespace

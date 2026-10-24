#include "EulerDomainStateRegistry.h"

#include <gtest/gtest.h>

#include <stdexcept>

namespace
{

using namespace ONEFLOW;

class TestState final : public EulerDomainState
{
public:
    int generation = 0;
};

class TestNs3DBackendState final : public Ns3DBackendState
{
public:
    TestNs3DBackendState(
        AccelBackendKind backendValue, int deviceIdValue,
        int * destructionsValue )
        : backend( backendValue ), deviceId( deviceIdValue ),
          destructions( destructionsValue )
    {
    }

    ~TestNs3DBackendState() override
    {
        if ( destructions != nullptr ) ++( *destructions );
    }

    AccelBackendKind Backend() const noexcept override
    {
        return backend;
    }

    int DeviceId() const noexcept override
    {
        return deviceId;
    }

private:
    AccelBackendKind backend;
    int deviceId;
    int * destructions;
};

EulerDomainStateKey Key(
    int solverIndex, int zone, int level, AccelBackendKind backend,
    int deviceId = -1 )
{
    return { solverIndex, zone, level, backend, deviceId };
}

TEST( EulerDomainStateRegistry, SeparatesSolverZoneGridBackendAndDevice )
{
    EulerDomainStateRegistry registry;
    const EulerDomainStateKey cpu = Key( 0, 2, 0, AccelBackendKind::CPU );
    const EulerDomainStateKey hip = Key( 0, 2, 0, AccelBackendKind::HIP, 0 );
    const EulerDomainStateKey otherDevice =
        Key( 0, 2, 0, AccelBackendKind::HIP, 1 );
    const EulerDomainStateKey otherZone = Key( 0, 3, 0, AccelBackendKind::CPU );

    registry.Insert( cpu, std::make_unique< TestState >() );
    registry.Insert( hip, std::make_unique< TestState >() );
    registry.Insert( otherDevice, std::make_unique< TestState >() );
    registry.Insert( otherZone, std::make_unique< TestState >() );

    EXPECT_EQ( registry.Size(), 4u );
    EXPECT_TRUE( registry.Contains( cpu ) );
    EXPECT_TRUE( registry.Contains( hip ) );
    EXPECT_TRUE( registry.Contains( otherDevice ) );
    EXPECT_TRUE( registry.Contains( otherZone ) );
    EXPECT_THROW(
        registry.Insert( cpu, std::make_unique< TestState >() ),
        std::logic_error );
}

TEST( EulerDomainStateRegistry, RejectsNullAndMissingStates )
{
    EulerDomainStateRegistry registry;
    const EulerDomainStateKey key = Key( 1, 1, 1, AccelBackendKind::CPU );

    EXPECT_THROW( registry.Insert( key, nullptr ), std::invalid_argument );
    EXPECT_THROW( registry.Get( key ), std::out_of_range );
    EXPECT_NO_THROW( registry.Erase( key ) );
}

TEST( EulerDomainStateRegistry, GetOrCreateReusesAndInvalidateReleases )
{
    EulerDomainStateRegistry registry;
    const EulerDomainStateKey key = Key( 3, 2, 1, AccelBackendKind::CPU );
    int creations = 0;

    EulerDomainState& first = registry.GetOrCreate( key, [&]() {
        ++creations;
        return std::make_unique<TestState>();
    } );
    EulerDomainState& second = registry.GetOrCreate( key, [&]() {
        ++creations;
        return std::make_unique<TestState>();
    } );

    EXPECT_EQ( &first, &second );
    EXPECT_EQ( creations, 1 );
    EXPECT_TRUE( registry.Invalidate( key ) );
    EXPECT_FALSE( registry.Contains( key ) );
    EXPECT_FALSE( registry.Invalidate( key ) );
}

TEST( EulerDomainStateRegistry, EraseAndClearReleaseOwnership )
{
    EulerDomainStateRegistry registry;
    const EulerDomainStateKey first = Key( 0, 0, 0, AccelBackendKind::CPU );
    const EulerDomainStateKey second = Key( 0, 1, 0, AccelBackendKind::CPU );
    registry.Insert( first, std::make_unique< TestState >() );
    registry.Insert( second, std::make_unique< TestState >() );

    registry.Erase( first );
    EXPECT_FALSE( registry.Contains( first ) );
    EXPECT_EQ( registry.Size(), 1u );
    registry.Clear();
    EXPECT_EQ( registry.Size(), 0u );
}


TEST( EulerDomainStateRegistry, RestartInvalidateCreatesFreshState )
{
    EulerDomainStateRegistry registry;
    const EulerDomainStateKey key = Key( 4, 1, 0, AccelBackendKind::CPU );
    int creations = 0;
    const auto factory = [&]() {
        auto state = std::make_unique< TestState >();
        state->generation = ++creations;
        return state;
    };

    EulerDomainState & first = registry.GetOrCreate( key, factory );
    EXPECT_EQ( static_cast< TestState & >( first ).generation, 1 );
    EXPECT_TRUE( registry.Invalidate( key ) );

    EulerDomainState & afterRestart = registry.GetOrCreate( key, factory );
    EXPECT_EQ( static_cast< TestState & >( afterRestart ).generation, 2 );
    EXPECT_EQ( creations, 2 );
}

TEST( Ns3DDeviceState, OwnsPersistentVerticalSliceBuffers )
{
    const EulerDomainProblem problem{ 4, 2, 5, 1.4, 0.01, 1.0, EulerDomainBoundary::Periodic };
    const EulerDomainStateKey key{ 2, 7, 0, AccelBackendKind::HIP, 0 };
    Ns3DDeviceState state( problem, key );
    state.ReserveFaces( 3 );
    std::vector< Real > values( 20, 1.0 );
    const EulerDomainConstFieldView field{ 4, 5, values.data() };
    state.Upload( field );

    EXPECT_EQ( state.key, key );
    EXPECT_EQ( state.nFaces, 3 );
    EXPECT_TRUE( state.uploaded );
    EXPECT_EQ( state.conservedState.size(), 20u );
    EXPECT_EQ( state.oldState.size(), 20u );
    EXPECT_EQ( state.rkScratch.size(), 20u );
    EXPECT_EQ( state.gradient.size(), 60u );
    EXPECT_EQ( state.faceFlux.size(), 15u );
    EXPECT_EQ( state.leftCell.size(), 3u );
    EXPECT_EQ( state.haloMetadata.size(), 3u );
}

TEST( Ns3DDeviceState, OwnsOpaqueBackendStateAndGenerationTokens )
{
    const EulerDomainProblem problem{
        4, 2, 5, 1.4, 0.01, 1.0, EulerDomainBoundary::Periodic };
    const EulerDomainStateKey key{ 2, 7, 0, AccelBackendKind::HIP, 0 };
    Ns3DDeviceState state( problem, key );
    const std::uint64_t initialTopology = state.TopologyGeneration();
    const std::uint64_t initialField = state.FieldGeneration();
    int destructions = 0;

    state.ReserveFaces( 3 );
    EXPECT_GT( state.TopologyGeneration(), initialTopology );
    state.AttachBackendState( std::make_unique< TestNs3DBackendState >(
        AccelBackendKind::HIP, 0, &destructions ) );
    EXPECT_TRUE( state.HasBackendState() );
    EXPECT_NE( state.OwnerToken(), nullptr );

    std::vector< Real > values( 20, 1.0 );
    state.Upload( EulerDomainConstFieldView{ 4, 5, values.data() } );
    EXPECT_GT( state.FieldGeneration(), initialField );

    state.ReleaseBackendState();
    EXPECT_FALSE( state.HasBackendState() );
    EXPECT_EQ( destructions, 1 );
    EXPECT_THROW(
        state.AttachBackendState( std::make_unique< TestNs3DBackendState >(
            AccelBackendKind::CPU, 0, &destructions ) ),
        std::invalid_argument );
    EXPECT_EQ( destructions, 2 );
    EXPECT_THROW(
        state.AttachBackendState( std::make_unique< TestNs3DBackendState >(
            AccelBackendKind::HIP, 1, &destructions ) ),
        std::invalid_argument );
    EXPECT_EQ( destructions, 3 );
}

TEST( Ns3DDeviceState, RegistryInvalidationReleasesBackendOwnership )
{
    const EulerDomainProblem problem{
        4, 2, 5, 1.4, 0.01, 1.0, EulerDomainBoundary::Periodic };
    const EulerDomainStateKey key{ 2, 7, 0, AccelBackendKind::HIP, 0 };
    EulerDomainStateRegistry registry;
    int destructions = 0;

    auto first = std::make_unique< Ns3DDeviceState >( problem, key );
    first->AttachBackendState( std::make_unique< TestNs3DBackendState >(
        AccelBackendKind::HIP, 0, &destructions ) );
    registry.Insert( key, std::move( first ) );
    EXPECT_EQ( destructions, 0 );
    EXPECT_TRUE( registry.Invalidate( key ) );
    EXPECT_EQ( destructions, 1 );

    auto restarted = std::make_unique< Ns3DDeviceState >( problem, key );
    restarted->AttachBackendState(
        std::make_unique< TestNs3DBackendState >(
            AccelBackendKind::HIP, 0, &destructions ) );
    registry.Insert( key, std::move( restarted ) );
    registry.Clear();
    EXPECT_EQ( destructions, 2 );
}

} // namespace

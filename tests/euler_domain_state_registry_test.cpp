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

EulerDomainStateKey Key(
    int solverIndex, int zone, int level, AccelBackendKind backend )
{
    return { solverIndex, zone, level, backend };
}

TEST( EulerDomainStateRegistry, SeparatesSolverZoneGridAndBackend )
{
    EulerDomainStateRegistry registry;
    const EulerDomainStateKey cpu = Key( 0, 2, 0, AccelBackendKind::CPU );
    const EulerDomainStateKey hip = Key( 0, 2, 0, AccelBackendKind::HIP );
    const EulerDomainStateKey otherZone = Key( 0, 3, 0, AccelBackendKind::CPU );

    registry.Insert( cpu, std::make_unique< TestState >() );
    registry.Insert( hip, std::make_unique< TestState >() );
    registry.Insert( otherZone, std::make_unique< TestState >() );

    EXPECT_EQ( registry.Size(), 3u );
    EXPECT_TRUE( registry.Contains( cpu ) );
    EXPECT_TRUE( registry.Contains( hip ) );
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

} // namespace

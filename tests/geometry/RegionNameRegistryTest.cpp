// tests/geometry/RegionNameRegistryTest.cpp
// RegionNameRegistry is instance-owned so tests do not touch process-wide
// RegionNameMap::Default() state (except when explicitly testing the facade).

#include <gtest/gtest.h>

#include "Boundary.h"
#include "Constant.h"

using namespace ONEFLOW;

TEST( RegionNameRegistryTest, AddAndFindAssignsStableIds )
{
    RegionNameRegistry registry;
    registry.AddRegion( "inlet" );
    registry.AddRegion( "outlet" );

    EXPECT_EQ( registry.FindRegionId( "inlet" ), 0 );
    EXPECT_EQ( registry.FindRegionId( "outlet" ), 1 );
    EXPECT_EQ( registry.FindRegionId( "missing" ), INVALID_INDEX );
}

TEST( RegionNameRegistryTest, DuplicateAddDoesNotChangeId )
{
    RegionNameRegistry registry;
    registry.AddRegion( "wall" );
    registry.AddRegion( "wall" );

    EXPECT_EQ( registry.FindRegionId( "wall" ), 0 );
    EXPECT_EQ( registry.GetNameMap().GetNameMap().size(), static_cast< std::size_t >( 1 ) );
}

TEST( RegionNameRegistryTest, ResetClearsEntries )
{
    RegionNameRegistry registry;
    registry.AddRegion( "farfield" );
    registry.Reset();

    EXPECT_EQ( registry.FindRegionId( "farfield" ), INVALID_INDEX );
    EXPECT_TRUE( registry.GetNameMap().GetNameMap().empty() );
}

TEST( RegionNameRegistryTest, InstancesAreIsolated )
{
    RegionNameRegistry a;
    RegionNameRegistry b;
    a.AddRegion( "a_only" );
    b.AddRegion( "b_only" );

    EXPECT_EQ( a.FindRegionId( "a_only" ), 0 );
    EXPECT_EQ( a.FindRegionId( "b_only" ), INVALID_INDEX );
    EXPECT_EQ( b.FindRegionId( "b_only" ), 0 );
    EXPECT_EQ( b.FindRegionId( "a_only" ), INVALID_INDEX );
}

TEST( RegionNameRegistryTest, DefaultFacadeForwardsToSharedRegistry )
{
    RegionNameMap::Reset();
    RegionNameMap::AddRegion( "default_region" );
    EXPECT_EQ( RegionNameMap::FindRegionId( "default_region" ), 0 );
    EXPECT_EQ( RegionNameMap::Default().FindRegionId( "default_region" ), 0 );
    RegionNameMap::Reset();
    EXPECT_EQ( RegionNameMap::FindRegionId( "default_region" ), INVALID_INDEX );
}

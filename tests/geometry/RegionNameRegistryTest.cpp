// tests/geometry/RegionNameRegistryTest.cpp
// RegionNameRegistry is instance-owned so tests do not touch process-wide
// RegionNameMap::Default() state (except when explicitly testing the facade).

#include <gtest/gtest.h>
#include <vector>
#include <string>

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

TEST( VolumeNameMapTest, DefaultFacadeIsIsolatedFromRegionNameMap )
{
    RegionNameMap::Reset();
    VolumeNameMap::Reset();

    RegionNameMap::AddRegion( "shared_spelling" );
    VolumeNameMap::AddRegion( "shared_spelling" );

    EXPECT_EQ( RegionNameMap::FindRegionId( "shared_spelling" ), 0 );
    EXPECT_EQ( VolumeNameMap::FindRegionId( "shared_spelling" ), 0 );

    RegionNameMap::AddRegion( "bc_only" );
    EXPECT_EQ( RegionNameMap::FindRegionId( "bc_only" ), 1 );
    EXPECT_EQ( VolumeNameMap::FindRegionId( "bc_only" ), INVALID_INDEX );

    VolumeNameMap::Reset();
    RegionNameMap::Reset();
}

TEST( VolumeNameMapTest, ResetClearsDefaultRegistry )
{
    VolumeNameMap::Reset();
    VolumeNameMap::AddRegion( "fluid" );
    EXPECT_EQ( VolumeNameMap::FindRegionId( "fluid" ), 0 );
    VolumeNameMap::Reset();
    EXPECT_EQ( VolumeNameMap::FindRegionId( "fluid" ), INVALID_INDEX );
}

TEST( RegisterRegionNamesTest, AssignsIdsInOrderAndDedups )
{
    RegionNameRegistry registry;
    const std::vector< std::string > names{ "inlet", "outlet", "inlet", "wall" };

    const IntField ids = RegisterRegionNames( registry, names );

    ASSERT_EQ( ids.size(), static_cast< std::size_t >( 4 ) );
    EXPECT_EQ( ids[ 0 ], 0 );
    EXPECT_EQ( ids[ 1 ], 1 );
    EXPECT_EQ( ids[ 2 ], 0 ); // duplicate name keeps first id
    EXPECT_EQ( ids[ 3 ], 2 );
    EXPECT_EQ( registry.GetNameMap().GetNameMap().size(), static_cast< std::size_t >( 3 ) );
}

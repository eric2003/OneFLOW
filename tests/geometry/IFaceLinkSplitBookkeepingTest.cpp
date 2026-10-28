// tests/geometry/IFaceLinkSplitBookkeepingTest.cpp
// Locks zone-scoped split bookkeeping accessors on IFaceLink.

#include <gtest/gtest.h>
#include <stdexcept>

#include "IFaceLink.h"
#include "GridHandles.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( IFaceLinkSplitBookkeepingTest, EnsureSplitBookkeepingGrowsNChild )
{
    Grids grids;
    IFaceLink link( grids );
    link.EnsureInterfaceZoneCount( 2 );
    EXPECT_EQ( link.GetSplitBookkeepingZoneCount(), static_cast< std::size_t >( 0 ) );

    link.EnsureSplitBookkeepingMatchesZones();
    ASSERT_EQ( link.GetSplitBookkeepingZoneCount(), static_cast< std::size_t >( 2 ) );
}

TEST( IFaceLinkSplitBookkeepingTest, MutableL2gNewBoundsCheck )
{
    Grids grids;
    IFaceLink link( grids );
    link.EnsureInterfaceZoneCount( 1 );
    // [[nodiscard]] refs: discard intentionally when only testing throw/no-throw.
    EXPECT_NO_THROW( static_cast< void >( link.MutableL2gNewForZone( 0 ) ) );
    EXPECT_THROW( static_cast< void >( link.MutableL2gNewForZone( 1 ) ), std::out_of_range );
    EXPECT_THROW( static_cast< void >( link.MutableL2gNewForZone( -1 ) ), std::out_of_range );
}

TEST( IFaceLinkSplitBookkeepingTest, GetL2gForZoneBoundsCheck )
{
    Grids grids;
    IFaceLink link( grids );
    link.EnsureInterfaceZoneCount( 1 );
    link.EnsureLocalInterfaceFaceCount( 0, 2 );
    // Slots exist but are unset (-1); register two local faces on new globals.
    link.RegisterLocalInterfaceFace( 0, true, 0, 0 );
    link.RegisterLocalInterfaceFace( 1, true, 0, 1 );
    EXPECT_EQ( link.GetL2gForZone( 0 ), ( IntField{ 0, 1 } ) );
    EXPECT_THROW( static_cast< void >( link.GetL2gForZone( 1 ) ), std::out_of_range );
}

TEST( IFaceLinkSplitBookkeepingTest, MutableNChildRequiresPriorEnsure )
{
    Grids grids;
    IFaceLink link( grids );
    link.EnsureInterfaceZoneCount( 1 );
    // nChild still empty -> out of range until Ensure
    EXPECT_THROW( static_cast< void >( link.MutableNChildForZone( 0 ) ), std::out_of_range );
    link.EnsureSplitBookkeepingMatchesZones();
    EXPECT_NO_THROW( static_cast< void >( link.MutableNChildForZone( 0 ) ) );
}


TEST( IFaceLinkSplitBookkeepingTest, EnsureInterfaceZoneCountSizesZoneTables )
{
    Grids grids;
    IFaceLink link( grids );
    EXPECT_TRUE( link.GetLocalToGlobalMap().empty() );
    link.EnsureInterfaceZoneCount( 3 );
    EXPECT_EQ( link.GetLocalToGlobalMap().size(), static_cast< std::size_t >( 3 ) );
    // Global-face tables unchanged
    EXPECT_TRUE( link.GetGlobalToZoneIds().empty() );
    EXPECT_TRUE( link.GetGlobalToLocalIds().empty() );
}

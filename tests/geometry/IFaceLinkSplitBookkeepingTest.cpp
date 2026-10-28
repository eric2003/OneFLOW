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
    link.l2g.resize( 2 );
    link.l2gNew.resize( 2 );
    EXPECT_TRUE( link.nChild.empty() );

    link.EnsureSplitBookkeepingMatchesZones();
    ASSERT_EQ( link.nChild.size(), static_cast< std::size_t >( 2 ) );
}

TEST( IFaceLinkSplitBookkeepingTest, MutableL2gNewBoundsCheck )
{
    Grids grids;
    IFaceLink link( grids );
    link.l2gNew.resize( 1 );
    EXPECT_NO_THROW( link.MutableL2gNewForZone( 0 ) );
    EXPECT_THROW( link.MutableL2gNewForZone( 1 ), std::out_of_range );
    EXPECT_THROW( link.MutableL2gNewForZone( -1 ), std::out_of_range );
}

TEST( IFaceLinkSplitBookkeepingTest, GetL2gForZoneBoundsCheck )
{
    Grids grids;
    IFaceLink link( grids );
    link.l2g.push_back( IntField{ 3, 4 } );
    EXPECT_EQ( link.GetL2gForZone( 0 ), ( IntField{ 3, 4 } ) );
    EXPECT_THROW( link.GetL2gForZone( 1 ), std::out_of_range );
}

TEST( IFaceLinkSplitBookkeepingTest, MutableNChildRequiresPriorEnsure )
{
    Grids grids;
    IFaceLink link( grids );
    link.l2g.resize( 1 );
    // nChild still empty -> out of range until Ensure
    EXPECT_THROW( link.MutableNChildForZone( 0 ), std::out_of_range );
    link.EnsureSplitBookkeepingMatchesZones();
    EXPECT_NO_THROW( link.MutableNChildForZone( 0 ) );
}


TEST( IFaceLinkSplitBookkeepingTest, EnsureInterfaceZoneCountSizesZoneTables )
{
    Grids grids;
    IFaceLink link( grids );
    EXPECT_TRUE( link.l2g.empty() );
    link.EnsureInterfaceZoneCount( 3 );
    EXPECT_EQ( link.l2g.size(), static_cast< std::size_t >( 3 ) );
    EXPECT_EQ( link.l2gNew.size(), static_cast< std::size_t >( 3 ) );
    // Global-face tables unchanged
    EXPECT_TRUE( link.gI2Zid.empty() );
    EXPECT_TRUE( link.g2l.empty() );
}

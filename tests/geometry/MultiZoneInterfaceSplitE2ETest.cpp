// tests/geometry/MultiZoneInterfaceSplitE2ETest.cpp
// Lightweight multi-zone interface split path without CGNS grids.
// Composes pure helpers already locked in isolation:
//   Register -> InitNew buffers -> RebuildZoneL2gNew -> BuildInverse (both zones)
//   -> ExpandBcRecord -> Commit.
// Intent: one readable scenario for the dual-map contract after a shared
// parent interface face splits into two children on both sides.

#include <gtest/gtest.h>

#include "FaceTopo.h"
#include "IFaceLink.h"
#include "GridHandles.h"
#include "HXDefine.h"
#include "Boundary.h"

using namespace ONEFLOW;

namespace {

void SeedTwoZoneTwoFaceMaps( IFaceLink & link )
{
    link.EnsureInterfaceZoneCount( 2 );
    link.EnsureLocalInterfaceFaceCount( 0, 2 );
    link.EnsureLocalInterfaceFaceCount( 1, 2 );
    link.RegisterLocalInterfaceFace( /*gIid=*/0, /*isNew=*/true,  /*zid=*/0, /*localFaceId=*/0 );
    link.RegisterLocalInterfaceFace( /*gIid=*/0, /*isNew=*/false, /*zid=*/1, /*localFaceId=*/0 );
    link.RegisterLocalInterfaceFace( /*gIid=*/1, /*isNew=*/true,  /*zid=*/0, /*localFaceId=*/1 );
    link.RegisterLocalInterfaceFace( /*gIid=*/1, /*isNew=*/false, /*zid=*/1, /*localFaceId=*/1 );
}

} // namespace

TEST( MultiZoneInterfaceSplitE2ETest, SplitOneSharedParentThenCommitBothZones )
{
    Grids grids;
    IFaceLink link( grids );
    SeedTwoZoneTwoFaceMaps( link );

    ASSERT_EQ( link.GetLocalToGlobalMap().size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( link.GetL2gForZone( 0 ), ( IntField{ 0, 1 } ) );
    EXPECT_EQ( link.GetL2gForZone( 1 ), ( IntField{ 0, 1 } ) );
    EXPECT_EQ( link.GetGlobalToZoneIds()[ 0 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( link.GetGlobalToLocalIds()[ 0 ], ( IntField{ 0, 0 } ) );

    link.InitNewLgMapping();
    link.EnsureSplitBookkeepingMatchesZones();

    const IntField parents{ 0, 1 };
    LinkField children( 2 );
    children[ 0 ] = IntField{ 2, 3 };
    children[ 1 ] = IntField{};

    for ( int zid = 0; zid < 2; ++ zid )
    {
        RebuildZoneL2gNewAfterInterfaceSplit(
            parents,
            children,
            link.MutableL2gNewForZone( zid ),
            link.MutableNChildForZone( zid ) );
    }

    EXPECT_EQ( link.MutableL2gNewForZone( 0 ), ( IntField{ 2, 3, 1 } ) );
    EXPECT_EQ( link.MutableL2gNewForZone( 1 ), ( IntField{ 2, 3, 1 } ) );
    EXPECT_EQ( link.MutableNChildForZone( 0 ), ( IntField{ 2, 0, 0, 0 } ) );

    BuildInverseInterfaceMapForZone(
        0, link.MutableL2gNewForZone( 0 ),
        link.MutableGlobalToZoneIdsNew(),
        link.MutableGlobalToLocalIdsNew() );
    BuildInverseInterfaceMapForZone(
        1, link.MutableL2gNewForZone( 1 ),
        link.MutableGlobalToZoneIdsNew(),
        link.MutableGlobalToLocalIdsNew() );

    ASSERT_GE( link.GetGlobalToZoneIdsNew().size(), static_cast< std::size_t >( 4 ) );
    EXPECT_EQ( link.GetGlobalToZoneIdsNew()[ 2 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( link.GetGlobalToLocalIdsNew()[ 2 ], ( IntField{ 0, 0 } ) );
    EXPECT_EQ( link.GetGlobalToZoneIdsNew()[ 3 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( link.GetGlobalToLocalIdsNew()[ 3 ], ( IntField{ 1, 1 } ) );
    EXPECT_EQ( link.GetGlobalToZoneIdsNew()[ 1 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( link.GetGlobalToLocalIdsNew()[ 1 ], ( IntField{ 2, 2 } ) );

    const IntField bcType{ BC::INTERFACE, BC::INTERFACE };
    const IntField bcNameId{ 100, 200 };
    const IntField childCounts{ 2, 0 };
    IntField outType, outName;
    ExpandBcRecordForSplitInterfaces( bcType, bcNameId, childCounts, outType, outName );
    EXPECT_EQ( outType, ( IntField{ BC::INTERFACE, BC::INTERFACE, BC::INTERFACE } ) );
    EXPECT_EQ( outName, ( IntField{ 100, 100, 200 } ) );

    link.UpdateLgMapping();

    EXPECT_EQ( link.GetL2gForZone( 0 ), ( IntField{ 2, 3, 1 } ) );
    EXPECT_EQ( link.GetL2gForZone( 1 ), ( IntField{ 2, 3, 1 } ) );
    EXPECT_EQ( link.GetGlobalToZoneIds()[ 2 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( link.GetGlobalToLocalIds()[ 2 ], ( IntField{ 0, 0 } ) );
    EXPECT_EQ( link.GetGlobalToZoneIds()[ 3 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( link.GetGlobalToLocalIds()[ 3 ], ( IntField{ 1, 1 } ) );

    // Legacy InitNew copies inverse rows: retired parent global 0 is not
    // stripped on Commit. No l2g entry points at it after the split.
    ASSERT_GE( link.GetGlobalToZoneIds().size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( link.GetGlobalToZoneIds()[ 0 ], ( IntField{ 0, 1 } ) );
}

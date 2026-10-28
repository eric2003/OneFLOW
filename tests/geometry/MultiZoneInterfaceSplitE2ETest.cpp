// tests/geometry/MultiZoneInterfaceSplitE2ETest.cpp
// Lightweight multi-zone interface split path without CGNS grids.
// Composes pure helpers already locked in isolation:
//   Register -> InitNew buffers -> RebuildZoneL2gNew -> BuildInverse (both zones)
//   -> ExpandBcRecord -> Commit.
// Intent: one readable scenario that states the dual-map contract after a
// single parent interface face splits into two children on both sides.

#include <gtest/gtest.h>

#include "FaceTopo.h"
#include "IFaceLink.h"
#include "GridHandles.h"
#include "HXDefine.h"
#include "Boundary.h"

using namespace ONEFLOW;

namespace {

// Two zones, two shared global interface faces (0 and 1). Face 0 will split
// into children 2 and 3; face 1 stays unsplit.
void SeedTwoZoneTwoFaceMaps( IFaceLink & link )
{
    link.EnsureInterfaceZoneCount( 2 );
    link.EnsureLocalInterfaceFaceCount( 0, 2 );
    link.EnsureLocalInterfaceFaceCount( 1, 2 );
    // global 0: zone0 local0 + zone1 local0
    link.RegisterLocalInterfaceFace( /*gIid=*/0, /*isNew=*/true,  /*zid=*/0, /*localFaceId=*/0 );
    link.RegisterLocalInterfaceFace( /*gIid=*/0, /*isNew=*/false, /*zid=*/1, /*localFaceId=*/0 );
    // global 1: zone0 local1 + zone1 local1
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

    // Stage *New tables from the live dual maps.
    link.InitNewLgMapping();
    link.EnsureSplitBookkeepingMatchesZones();

    // Both zones see the same parent list and the same child replacement for
    // parent global face 0. Parent global face 1 stays as-is.
    const IntField parents{ 0, 1 };
    LinkField children( 2 );
    children[ 0 ] = IntField{ 2, 3 }; // split
    children[ 1 ] = IntField{};      // keep

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
    // parent0 count 2 + trailing 0 per child + parent1 count 0
    EXPECT_EQ( link.MutableNChildForZone( 0 ), ( IntField{ 2, 0, 0, 0 } ) );

    // Inverse maps: both zones publish local indices for globals 2, 3, 1.
    BuildInverseInterfaceMapForZone(
        0, link.MutableL2gNewForZone( 0 ),
        link.MutableGlobalToZoneIdsNew(),
        link.MutableGlobalToLocalIdsNew() );
    BuildInverseInterfaceMapForZone(
        1, link.MutableL2gNewForZone( 1 ),
        link.MutableGlobalToZoneIdsNew(),
        link.MutableGlobalToLocalIdsNew() );

    // Global 2 and 3 must list both zones; local slots are 0 and 1 after rebuild.
    ASSERT_GE( link.GetGlobalToZoneIdsNew().size(), static_cast< std::size_t >( 4 ) );
    EXPECT_EQ( link.GetGlobalToZoneIdsNew()[ 2 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( link.GetGlobalToLocalIdsNew()[ 2 ], ( IntField{ 0, 0 } ) );
    EXPECT_EQ( link.GetGlobalToZoneIdsNew()[ 3 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( link.GetGlobalToLocalIdsNew()[ 3 ], ( IntField{ 1, 1 } ) );
    // Unsplit parent global 1 keeps two participants at local index 2.
    EXPECT_EQ( link.GetGlobalToZoneIdsNew()[ 1 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( link.GetGlobalToLocalIdsNew()[ 1 ], ( IntField{ 2, 2 } ) );

    // BC expand is zone-local: one interface face splits into two BC slots.
    const IntField bcType{ BC::INTERFACE, BC::INTERFACE };
    const IntField bcNameId{ 100, 200 };
    const IntField childCounts{ 2, 0 }; // parallel to interface-local order
    IntField outType, outName;
    ExpandBcRecordForSplitInterfaces( bcType, bcNameId, childCounts, outType, outName );
    EXPECT_EQ( outType, ( IntField{ BC::INTERFACE, BC::INTERFACE, BC::INTERFACE } ) );
    EXPECT_EQ( outName, ( IntField{ 100, 100, 200 } ) );

    // Commit staged tables into the live dual maps.
    link.UpdateLgMapping();

    EXPECT_EQ( link.GetL2gForZone( 0 ), ( IntField{ 2, 3, 1 } ) );
    EXPECT_EQ( link.GetL2gForZone( 1 ), ( IntField{ 2, 3, 1 } ) );
    EXPECT_EQ( link.GetGlobalToZoneIds()[ 2 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( link.GetGlobalToLocalIds()[ 2 ], ( IntField{ 0, 0 } ) );
    EXPECT_EQ( link.GetGlobalToZoneIds()[ 3 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( link.GetGlobalToLocalIds()[ 3 ], ( IntField{ 1, 1 } ) );

    // Legacy InitNew copies inverse rows: the retired parent global 0 is not
    // stripped on Commit. No l2g entry points at it after the split.
    ASSERT_GE( link.GetGlobalToZoneIds().size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( link.GetGlobalToZoneIds()[ 0 ], ( IntField{ 0, 1 } ) );
}

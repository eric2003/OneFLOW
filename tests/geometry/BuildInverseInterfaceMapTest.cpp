// tests/geometry/BuildInverseInterfaceMapTest.cpp
// Locks inverse map build: l2gNew[local] -> gI2ZidNew/g2lNew[global].

#include <gtest/gtest.h>
#include <stdexcept>

#include "FaceTopo.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( BuildInverseInterfaceMapTest, RecordsZoneAndLocalOnGlobalFace )
{
    LinkField gI2ZidNew;
    LinkField g2lNew;
    const IntField l2gNewForZone{ 0, 2 }; // local 0 -> global 0, local 1 -> global 2

    BuildInverseInterfaceMapForZone( /*zoneId=*/3, l2gNewForZone, gI2ZidNew, g2lNew );

    ASSERT_EQ( gI2ZidNew.size(), static_cast< std::size_t >( 3 ) );
    ASSERT_EQ( g2lNew.size(), static_cast< std::size_t >( 3 ) );
    EXPECT_EQ( gI2ZidNew[ 0 ], ( IntField{ 3 } ) );
    EXPECT_EQ( g2lNew[ 0 ], ( IntField{ 0 } ) );
    EXPECT_TRUE( gI2ZidNew[ 1 ].empty() );
    EXPECT_EQ( gI2ZidNew[ 2 ], ( IntField{ 3 } ) );
    EXPECT_EQ( g2lNew[ 2 ], ( IntField{ 1 } ) );
}

TEST( BuildInverseInterfaceMapTest, SecondZoneFillsPairSlot )
{
    LinkField gI2ZidNew;
    LinkField g2lNew;
    BuildInverseInterfaceMapForZone( 0, IntField{ 0 }, gI2ZidNew, g2lNew );
    BuildInverseInterfaceMapForZone( 1, IntField{ 0 }, gI2ZidNew, g2lNew );

    ASSERT_EQ( gI2ZidNew[ 0 ].size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( gI2ZidNew[ 0 ][ 0 ], 0 );
    EXPECT_EQ( gI2ZidNew[ 0 ][ 1 ], 1 );
    EXPECT_EQ( g2lNew[ 0 ][ 0 ], 0 );
    EXPECT_EQ( g2lNew[ 0 ][ 1 ], 0 );
}

TEST( BuildInverseInterfaceMapTest, RefreshesExistingZoneSlotDuringRebuild )
{
    LinkField gI2ZidNew{ IntField{ 0, 1 } };
    LinkField g2lNew{ IntField{ 4, 8 } };

    // Global face 0 moved to local slot 1 in zone 0 after rebuilding.
    BuildInverseInterfaceMapForZone( 0, IntField{ 5, 0 }, gI2ZidNew, g2lNew );

    EXPECT_EQ( gI2ZidNew[ 0 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( g2lNew[ 0 ], ( IntField{ 1, 8 } ) );
    EXPECT_EQ( gI2ZidNew[ 5 ], ( IntField{ 0 } ) );
    EXPECT_EQ( g2lNew[ 5 ], ( IntField{ 0 } ) );
}

TEST( BuildInverseInterfaceMapTest, DoesNotDuplicateExistingZoneParticipant )
{
    LinkField gI2ZidNew{ IntField{ 0 } };
    LinkField g2lNew{ IntField{ 6 } };

    BuildInverseInterfaceMapForZone( 0, IntField{ 0 }, gI2ZidNew, g2lNew );

    EXPECT_EQ( gI2ZidNew[ 0 ], ( IntField{ 0 } ) );
    EXPECT_EQ( g2lNew[ 0 ], ( IntField{ 0 } ) );
}

TEST( BuildInverseInterfaceMapTest, LegacyPolicyIgnoresThirdZoneParticipant )
{
    LinkField gI2ZidNew;
    LinkField g2lNew;
    BuildInverseInterfaceMapForZone( 0, IntField{ 0 }, gI2ZidNew, g2lNew );
    BuildInverseInterfaceMapForZone( 1, IntField{ 0 }, gI2ZidNew, g2lNew );

    const LinkField zoneIdsBeforeThirdZone = gI2ZidNew;
    const LinkField localIdsBeforeThirdZone = g2lNew;
    BuildInverseInterfaceMapForZone( 2, IntField{ 0 }, gI2ZidNew, g2lNew );

    EXPECT_EQ( gI2ZidNew, zoneIdsBeforeThirdZone );
    EXPECT_EQ( g2lNew, localIdsBeforeThirdZone );
}

TEST( BuildInverseInterfaceMapTest, RejectsNegativeZoneIdWithoutChangingMaps )
{
    LinkField gI2ZidNew{ IntField{ 4 } };
    LinkField g2lNew{ IntField{ 7 } };

    EXPECT_THROW(
        BuildInverseInterfaceMapForZone( -1, IntField{ 0 }, gI2ZidNew, g2lNew ),
        std::out_of_range );

    EXPECT_EQ( gI2ZidNew, ( LinkField{ IntField{ 4 } } ) );
    EXPECT_EQ( g2lNew, ( LinkField{ IntField{ 7 } } ) );
}

TEST( BuildInverseInterfaceMapTest, RejectsNegativeGlobalIdWithoutPartialMutation )
{
    LinkField gI2ZidNew{ IntField{ 8 } };
    LinkField g2lNew{ IntField{ 9 } };
    const LinkField originalZoneIds = gI2ZidNew;
    const LinkField originalLocalIds = g2lNew;

    // A valid, out-of-range ID comes first: a single-pass implementation would
    // resize and append before discovering the invalid second ID.
    EXPECT_THROW(
        BuildInverseInterfaceMapForZone(
            0, IntField{ 2, -1 }, gI2ZidNew, g2lNew ),
        std::out_of_range );

    EXPECT_EQ( gI2ZidNew, originalZoneIds );
    EXPECT_EQ( g2lNew, originalLocalIds );
}

TEST( BuildInverseInterfaceMapTest, RejectsMismatchedTopLevelTables )
{
    LinkField gI2ZidNew;
    LinkField g2lNew;
    gI2ZidNew.resize( 1 );
    // g2lNew empty -> size mismatch
    EXPECT_THROW(
        BuildInverseInterfaceMapForZone( 0, IntField{ 0 }, gI2ZidNew, g2lNew ),
        std::logic_error );
}

TEST( BuildInverseInterfaceMapTest, RejectsDuplicateGlobalFaceWithinZoneWithoutMutation )
{
    LinkField gI2ZidNew{ IntField{ 4 } };
    LinkField g2lNew{ IntField{ 7 } };
    const LinkField originalZoneIds = gI2ZidNew;
    const LinkField originalLocalIds = g2lNew;

    EXPECT_THROW(
        BuildInverseInterfaceMapForZone(
            2, IntField{ 1, 1 }, gI2ZidNew, g2lNew ),
        std::logic_error );

    EXPECT_EQ( gI2ZidNew, originalZoneIds );
    EXPECT_EQ( g2lNew, originalLocalIds );
}

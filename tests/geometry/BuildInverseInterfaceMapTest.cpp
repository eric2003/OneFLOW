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

TEST( BuildInverseInterfaceMapTest, DoesNotExceedTwoZoneSlots )
{
    LinkField gI2ZidNew;
    LinkField g2lNew;
    BuildInverseInterfaceMapForZone( 0, IntField{ 0 }, gI2ZidNew, g2lNew );
    BuildInverseInterfaceMapForZone( 1, IntField{ 0 }, gI2ZidNew, g2lNew );
    BuildInverseInterfaceMapForZone( 2, IntField{ 0 }, gI2ZidNew, g2lNew );

    EXPECT_EQ( gI2ZidNew[ 0 ].size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( gI2ZidNew[ 0 ][ 0 ], 0 );
    EXPECT_EQ( gI2ZidNew[ 0 ][ 1 ], 1 );
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

TEST( BuildInverseInterfaceMapTest, RejectsNegativeGlobalId )
{
    LinkField gI2ZidNew;
    LinkField g2lNew;
    EXPECT_THROW(
        BuildInverseInterfaceMapForZone( 0, IntField{ -1 }, gI2ZidNew, g2lNew ),
        std::out_of_range );
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

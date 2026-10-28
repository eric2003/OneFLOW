// tests/geometry/RegisterLocalInterfaceFaceTest.cpp
// Locks gI2Zid / g2l / l2g registration without constructing IFaceLink.

#include <gtest/gtest.h>
#include <stdexcept>

#include "IFaceLink.h"
#include "HXDefine.h"

using namespace ONEFLOW;

namespace {

void SeedTwoZonesThreeLocalSlots( LinkField & l2g )
{
    l2g.clear();
    l2g.push_back( IntField{ -1, -1, -1 } ); // zone 0
    l2g.push_back( IntField{ -1, -1, -1 } ); // zone 1
}

} // namespace

TEST( RegisterLocalInterfaceFaceTest, NewGlobalFaceAppendsParallelRows )
{
    LinkField l2g;
    SeedTwoZonesThreeLocalSlots( l2g );
    LinkField gI2Zid;
    LinkField g2l;

    RegisterLocalInterfaceFaceOnMaps(
        /*gIid=*/0, /*isNew=*/true, /*zid=*/0, /*localFaceId=*/1,
        l2g, gI2Zid, g2l );

    EXPECT_EQ( l2g[ 0 ][ 1 ], 0 );
    ASSERT_EQ( gI2Zid.size(), static_cast< std::size_t >( 1 ) );
    ASSERT_EQ( g2l.size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( gI2Zid[ 0 ], ( IntField{ 0 } ) );
    EXPECT_EQ( g2l[ 0 ], ( IntField{ 1 } ) );
}

TEST( RegisterLocalInterfaceFaceTest, ExistingGlobalFacePushesSecondSide )
{
    LinkField l2g;
    SeedTwoZonesThreeLocalSlots( l2g );
    LinkField gI2Zid;
    LinkField g2l;

    RegisterLocalInterfaceFaceOnMaps( 0, true, 0, 0, l2g, gI2Zid, g2l );
    RegisterLocalInterfaceFaceOnMaps( 0, false, 1, 2, l2g, gI2Zid, g2l );

    EXPECT_EQ( l2g[ 0 ][ 0 ], 0 );
    EXPECT_EQ( l2g[ 1 ][ 2 ], 0 );
    ASSERT_EQ( gI2Zid[ 0 ].size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( gI2Zid[ 0 ][ 0 ], 0 );
    EXPECT_EQ( gI2Zid[ 0 ][ 1 ], 1 );
    EXPECT_EQ( g2l[ 0 ][ 0 ], 0 );
    EXPECT_EQ( g2l[ 0 ][ 1 ], 2 );
}


TEST( RegisterLocalInterfaceFaceTest, RejectsThirdZoneParticipantWithoutChangingMaps )
{
    LinkField l2g;
    SeedTwoZonesThreeLocalSlots( l2g );
    l2g.push_back( IntField{ -1, -1, -1 } ); // zone 2
    LinkField gI2Zid;
    LinkField g2l;

    RegisterLocalInterfaceFaceOnMaps( 0, true, 0, 0, l2g, gI2Zid, g2l );
    RegisterLocalInterfaceFaceOnMaps( 0, false, 1, 0, l2g, gI2Zid, g2l );
    const LinkField originalL2g = l2g;
    const LinkField originalZoneIds = gI2Zid;
    const LinkField originalLocalIds = g2l;

    EXPECT_THROW(
        RegisterLocalInterfaceFaceOnMaps( 0, false, 2, 0, l2g, gI2Zid, g2l ),
        std::logic_error );

    EXPECT_EQ( l2g, originalL2g );
    EXPECT_EQ( gI2Zid, originalZoneIds );
    EXPECT_EQ( g2l, originalLocalIds );
}

TEST( RegisterLocalInterfaceFaceTest, RejectsDuplicateLocalParticipantWithoutChangingMaps )
{
    LinkField l2g;
    SeedTwoZonesThreeLocalSlots( l2g );
    LinkField gI2Zid;
    LinkField g2l;

    RegisterLocalInterfaceFaceOnMaps( 0, true, 0, 1, l2g, gI2Zid, g2l );
    const LinkField originalL2g = l2g;
    const LinkField originalZoneIds = gI2Zid;
    const LinkField originalLocalIds = g2l;

    EXPECT_THROW(
        RegisterLocalInterfaceFaceOnMaps( 0, false, 0, 1, l2g, gI2Zid, g2l ),
        std::logic_error );

    EXPECT_EQ( l2g, originalL2g );
    EXPECT_EQ( gI2Zid, originalZoneIds );
    EXPECT_EQ( g2l, originalLocalIds );
}

TEST( RegisterLocalInterfaceFaceTest, RejectsBadNewGlobalId )
{
    LinkField l2g;
    SeedTwoZonesThreeLocalSlots( l2g );
    LinkField gI2Zid;
    LinkField g2l;

    EXPECT_THROW(
        RegisterLocalInterfaceFaceOnMaps( 1, true, 0, 0, l2g, gI2Zid, g2l ),
        std::logic_error );
}

TEST( RegisterLocalInterfaceFaceTest, RejectsOutOfRangeLocalSlot )
{
    LinkField l2g;
    SeedTwoZonesThreeLocalSlots( l2g );
    LinkField gI2Zid;
    LinkField g2l;

    EXPECT_THROW(
        RegisterLocalInterfaceFaceOnMaps( 0, true, 0, 9, l2g, gI2Zid, g2l ),
        std::out_of_range );
}

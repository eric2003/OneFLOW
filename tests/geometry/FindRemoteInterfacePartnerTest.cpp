// tests/geometry/FindRemoteInterfacePartnerTest.cpp
// Locks remote zone/local selection on a global interface face.

#include <gtest/gtest.h>
#include <stdexcept>

#include "IFaceLink.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( FindRemoteInterfacePartnerTest, PicksNonSelfEntry )
{
    LinkField zoneIds;
    LinkField localIds;
    zoneIds.push_back( IntField{ 0, 1 } );
    localIds.push_back( IntField{ 3, 7 } );

    int rz = -1, rl = -1;
    EXPECT_TRUE(
        FindRemoteInterfacePartner( 0, /*selfZone=*/0, /*selfLocal=*/3,
                                    zoneIds, localIds, rz, rl ) );
    EXPECT_EQ( rz, 1 );
    EXPECT_EQ( rl, 7 );
}

TEST( FindRemoteInterfacePartnerTest, SelfOnSecondSlotStillFindsFirst )
{
    LinkField zoneIds;
    LinkField localIds;
    zoneIds.push_back( IntField{ 2, 4 } );
    localIds.push_back( IntField{ 1, 9 } );

    int rz = -1, rl = -1;
    EXPECT_TRUE(
        FindRemoteInterfacePartner( 0, 4, 9, zoneIds, localIds, rz, rl ) );
    EXPECT_EQ( rz, 2 );
    EXPECT_EQ( rl, 1 );
}

TEST( FindRemoteInterfacePartnerTest, ReturnsFalseWhenOnlySelfPresent )
{
    LinkField zoneIds;
    LinkField localIds;
    zoneIds.push_back( IntField{ 0 } );
    localIds.push_back( IntField{ 5 } );

    int rz = -1, rl = -1;
    EXPECT_FALSE(
        FindRemoteInterfacePartner( 0, 0, 5, zoneIds, localIds, rz, rl ) );
}

TEST( FindRemoteInterfacePartnerTest, RejectsBadGlobalId )
{
    LinkField zoneIds;
    LinkField localIds;
    zoneIds.push_back( IntField{ 0 } );
    localIds.push_back( IntField{ 0 } );
    int rz = -1, rl = -1;
    EXPECT_THROW(
        FindRemoteInterfacePartner( 1, 0, 0, zoneIds, localIds, rz, rl ),
        std::logic_error );
}

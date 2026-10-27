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

TEST( FindRemoteInterfacePartnerTest, DoesNotReturnUnrelatedParticipant )
{
    LinkField zoneIds;
    LinkField localIds;
    zoneIds.push_back( IntField{ 2, 4 } );
    localIds.push_back( IntField{ 1, 9 } );

    int remoteZone = 77;
    int remoteLocal = 88;
    EXPECT_FALSE(
        FindRemoteInterfacePartner(
            0, /*selfZone=*/0, /*selfLocal=*/3,
            zoneIds, localIds, remoteZone, remoteLocal ) );
    EXPECT_EQ( remoteZone, 77 );
    EXPECT_EQ( remoteLocal, 88 );
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

TEST( FindPeriodicInterfacePartnerTest, SkipsSelfAndSelectsFirstOtherParticipant )
{
    const IntField zoneIds = { 0, 2, 3 };
    const IntField localIds = { 4, 8, 9 };
    int remoteZone = -1;
    int remoteLocal = -1;

    EXPECT_TRUE( FindPeriodicInterfacePartner(
        0, 4, zoneIds, localIds, remoteZone, remoteLocal ) );
    EXPECT_EQ( remoteZone, 2 );
    EXPECT_EQ( remoteLocal, 8 );
}

TEST( FindPeriodicInterfacePartnerTest, ReturnsFalseForSelfOnlyAndPreservesOutputs )
{
    const IntField zoneIds = { 0 };
    const IntField localIds = { 4 };
    int remoteZone = 77;
    int remoteLocal = 88;

    EXPECT_FALSE( FindPeriodicInterfacePartner(
        0, 4, zoneIds, localIds, remoteZone, remoteLocal ) );
    EXPECT_EQ( remoteZone, 77 );
    EXPECT_EQ( remoteLocal, 88 );
}

TEST( FindPeriodicInterfacePartnerTest, RejectsMisalignedRows )
{
    const IntField zoneIds = { 0, 2 };
    const IntField localIds = { 4 };
    int remoteZone = 77;
    int remoteLocal = 88;

    EXPECT_THROW( FindPeriodicInterfacePartner(
        0, 4, zoneIds, localIds, remoteZone, remoteLocal ), std::logic_error );
    EXPECT_EQ( remoteZone, 77 );
    EXPECT_EQ( remoteLocal, 88 );
}

TEST( FindPeriodicInterfacePartnerTest, UsesFirstParticipantWhenCallerIsNotInPartnerRow )
{
    // Unlike topological lookup, geometric periodic matching searches a
    // different face row, so the caller may not be a participant in that row.
    const IntField zoneIds = { 5, 7 };
    const IntField localIds = { 12, 14 };
    int remoteZone = 77;
    int remoteLocal = 88;

    EXPECT_TRUE( FindPeriodicInterfacePartner(
        0, 4, zoneIds, localIds, remoteZone, remoteLocal ) );
    EXPECT_EQ( remoteZone, 5 );
    EXPECT_EQ( remoteLocal, 12 );
}

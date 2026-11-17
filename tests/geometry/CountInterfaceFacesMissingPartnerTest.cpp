// tests/geometry/CountInterfaceFacesMissingPartnerTest.cpp
// Locks unpaired global-face counting used by MatchInterfaceTopology.

#include <gtest/gtest.h>
#include <stdexcept>

#include "IFaceLink.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( CountInterfaceFacesMissingPartnerTest, CountsZeroAndOneZoneFaces )
{
    LinkField zoneIds;
    LinkField localIds;
    zoneIds.push_back( IntField{ 0 } );
    localIds.push_back( IntField{ 0 } );
    zoneIds.push_back( IntField{ 0, 1 } );
    localIds.push_back( IntField{ 1, 0 } );
    zoneIds.push_back( IntField{ 0, 1, 2 } );
    localIds.push_back( IntField{ 2, 1, 0 } );

    const IntField l2g{ 0, 1, 2 };
    EXPECT_EQ(
        CountInterfaceFacesMissingPartner( l2g, zoneIds, localIds ), 1 );
}

TEST( CountInterfaceFacesMissingPartnerTest, EmptyLocalListIsZero )
{
    LinkField zoneIds;
    LinkField localIds;
    EXPECT_EQ(
        CountInterfaceFacesMissingPartner( IntField{}, zoneIds, localIds ), 0 );
}

TEST( CountInterfaceFacesMissingPartnerTest, RejectsInconsistentTables )
{
    LinkField zoneIds;
    LinkField localIds;
    zoneIds.push_back( IntField{ 0, 1 } );
    localIds.push_back( IntField{ 0 } );
    EXPECT_THROW(
        CountInterfaceFacesMissingPartner( IntField{ 0 }, zoneIds, localIds ),
        std::logic_error );
}

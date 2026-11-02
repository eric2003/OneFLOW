// tests/geometry/GlobalFaceMappingTest.cpp
// Locks ValidateGlobalFaceMapping(s) without constructing IFaceLink or grids.

#include <gtest/gtest.h>
#include <stdexcept>

#include "IFaceLink.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( GlobalFaceMappingTest, AcceptsAlignedZoneAndLocalLists )
{
    LinkField zoneIds;
    LinkField localIds;
    zoneIds.push_back( IntField{ 0, 1 } );
    localIds.push_back( IntField{ 3, 4 } );

    EXPECT_NO_THROW(
        ValidateGlobalFaceMapping( 0, zoneIds, localIds, "test" ) );
    EXPECT_NO_THROW(
        ValidateGlobalFaceMappings( zoneIds, localIds, "test" ) );
}

TEST( GlobalFaceMappingTest, RejectsOutOfRangeGlobalId )
{
    LinkField zoneIds;
    LinkField localIds;
    zoneIds.push_back( IntField{ 0 } );
    localIds.push_back( IntField{ 1 } );

    EXPECT_THROW(
        ValidateGlobalFaceMapping( 1, zoneIds, localIds, "test" ),
        std::logic_error );
    EXPECT_THROW(
        ValidateGlobalFaceMapping( -1, zoneIds, localIds, "test" ),
        std::logic_error );
}

TEST( GlobalFaceMappingTest, RejectsMismatchedZoneLocalLengths )
{
    LinkField zoneIds;
    LinkField localIds;
    zoneIds.push_back( IntField{ 0, 1 } );
    localIds.push_back( IntField{ 3 } ); // different length

    EXPECT_THROW(
        ValidateGlobalFaceMapping( 0, zoneIds, localIds, "test" ),
        std::logic_error );
}

TEST( GlobalFaceMappingTest, RejectsTableSizeMismatch )
{
    LinkField zoneIds;
    LinkField localIds;
    zoneIds.push_back( IntField{ 0 } );
    zoneIds.push_back( IntField{ 1 } );
    localIds.push_back( IntField{ 0 } ); // only one entry

    EXPECT_THROW(
        ValidateGlobalFaceMappings( zoneIds, localIds, "test" ),
        std::logic_error );
}

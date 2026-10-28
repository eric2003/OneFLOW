// tests/geometry/RequireSoleZoneIdForGlobalFaceTest.cpp
// Locks line-split policy: global face must list exactly one zone.

#include <gtest/gtest.h>
#include <stdexcept>

#include "FaceSearch.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( RequireSoleZoneIdForGlobalFaceTest, ReturnsTheSoleZone )
{
    LinkField zones;
    zones.push_back( IntField{ 3 } );
    EXPECT_EQ( RequireSoleZoneIdForGlobalFace( zones, 0 ), 3 );
}

TEST( RequireSoleZoneIdForGlobalFaceTest, RejectsZeroOrManyZones )
{
    LinkField zones;
    zones.push_back( IntField{} );
    zones.push_back( IntField{ 0, 1 } );
    EXPECT_THROW( RequireSoleZoneIdForGlobalFace( zones, 0 ), std::logic_error );
    EXPECT_THROW( RequireSoleZoneIdForGlobalFace( zones, 1 ), std::logic_error );
}

TEST( RequireSoleZoneIdForGlobalFaceTest, RejectsBadGlobalId )
{
    LinkField zones;
    zones.push_back( IntField{ 0 } );
    EXPECT_THROW( RequireSoleZoneIdForGlobalFace( zones, 1 ), std::out_of_range );
    EXPECT_THROW( RequireSoleZoneIdForGlobalFace( zones, -1 ), std::out_of_range );
}

// tests/geometry/InterfaceChildFaceIdsTest.cpp
// GetInterfaceChildFaceIds (table overload) maps local interface face ->
// global parent -> child face id list, and validates child node tables.

#include <gtest/gtest.h>
#include <stdexcept>

#include "FaceTopo.h"
#include "FaceSearch.h"
#include "HXDefine.h"

using namespace ONEFLOW;

namespace {

// zone 0: local interface faces 0,1 map to global parents 0,1
// parent 0 has children {0,1}; parent 1 has no children
void SeedSplitTable( LinkField & l2g, InterfaceFaceSplitTable & split )
{
    l2g.clear();
    l2g.push_back( IntField{ 0, 1 } );

    split.Clear();
    split.childFaceIdsOfParent.resize( 2 );
    split.childFaceIdsOfParent[ 0 ] = IntField{ 0, 1 };
    split.childFaceIdsOfParent[ 1 ] = IntField{};

    // child faces 0 and 1: matching node id / flag lengths
    split.childFaceNodeIds.resize( 2 );
    split.childFaceNodeFlags.resize( 2 );
    split.childFaceNodeIds[ 0 ] = IntField{ 10, 11, 12 };
    split.childFaceNodeFlags[ 0 ] = IntField{ 0, 0, 1 };
    split.childFaceNodeIds[ 1 ] = IntField{ 11, 12, 13 };
    split.childFaceNodeFlags[ 1 ] = IntField{ 0, 1, 0 };
}

} // namespace

TEST( InterfaceChildFaceIdsTest, ReturnsChildrenForSplitParent )
{
    LinkField l2g;
    InterfaceFaceSplitTable split;
    SeedSplitTable( l2g, split );

    const IntField & children =
        GetInterfaceChildFaceIds( l2g, split, /*zoneId=*/0, /*localInterfaceFaceId=*/0 );

    ASSERT_EQ( children.size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( children[ 0 ], 0 );
    EXPECT_EQ( children[ 1 ], 1 );
}

TEST( InterfaceChildFaceIdsTest, UnsplitParentReturnsEmptyList )
{
    LinkField l2g;
    InterfaceFaceSplitTable split;
    SeedSplitTable( l2g, split );

    const IntField & children =
        GetInterfaceChildFaceIds( l2g, split, 0, /*localInterfaceFaceId=*/1 );

    EXPECT_TRUE( children.empty() );
}

TEST( InterfaceChildFaceIdsTest, RejectsOutOfRangeLocalFace )
{
    LinkField l2g;
    InterfaceFaceSplitTable split;
    SeedSplitTable( l2g, split );

    EXPECT_THROW(
        GetInterfaceChildFaceIds( l2g, split, 0, /*localInterfaceFaceId=*/5 ),
        std::out_of_range );
}

TEST( InterfaceChildFaceIdsTest, ChildIdentityLookupDoesNotValidateLegacyNodeMappings )
{
    LinkField l2g;
    InterfaceFaceSplitTable split;
    SeedSplitTable( l2g, split );
    // Node mapping validity belongs to the mapping collector, not ID lookup.
    split.childFaceNodeFlags[ 0 ] = IntField{ 0 };

    const IntField & children = GetInterfaceChildFaceIds( l2g, split, 0, 0 );

    EXPECT_EQ( children, ( IntField{ 0, 1 } ) );
}

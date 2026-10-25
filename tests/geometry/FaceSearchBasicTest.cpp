// tests/geometry/FaceSearchBasicTest.cpp
// FaceSearchBasic pairs undirected identity (HXLookup) with directed storage
// (faceArray keeps the node order of the first AddFace).

#include <gtest/gtest.h>
#include <stdexcept>

#include "FaceSearch.h"
#include "Constant.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( FaceSearchBasicTest, AddFaceAssignsSequentialIds )
{
    FaceSearchBasic search;
    EXPECT_EQ( search.AddFace( IntField{ 1, 2, 3 } ), 0 );
    EXPECT_EQ( search.AddFace( IntField{ 4, 5, 6 } ), 1 );
    EXPECT_EQ( search.Size(), static_cast< std::size_t >( 2 ) );
    ASSERT_EQ( search.Size(), static_cast< std::size_t >( 2 ) );
}

TEST( FaceSearchBasicTest, DuplicateNodesSameIdKeepsFirstOrder )
{
    // Intent: identity ignores order; orientation memory keeps first ring.
    FaceSearchBasic search;
    const int id0 = search.AddFace( IntField{ 1, 2, 3 } );
    const int id1 = search.AddFace( IntField{ 3, 1, 2 } );

    EXPECT_EQ( id0, id1 );
    EXPECT_EQ( search.Size(), static_cast< std::size_t >( 1 ) );
    ASSERT_EQ( search.Size(), static_cast< std::size_t >( 1 ) );
    EXPECT_EQ( search.GetFaceNodeRing( 0 ), ( IntField{ 1, 2, 3 } ) );
}

TEST( FaceSearchBasicTest, FindFaceIsOrderInvariant )
{
    FaceSearchBasic search;
    search.AddFace( IntField{ 2, 0, 1 } );
    EXPECT_EQ( search.FindFace( IntField{ 0, 1, 2 } ), 0 );
    EXPECT_EQ( search.FindFace( IntField{ 9, 8, 7 } ), INVALID_INDEX );
}

TEST( FaceSearchBasicTest, ClearResetsLookupAndFaceArray )
{
    FaceSearchBasic search;
    search.AddFace( IntField{ 1, 2 } );
    search.Clear();
    EXPECT_EQ( search.Size(), static_cast< std::size_t >( 0 ) );
    EXPECT_EQ( search.Size(), static_cast< std::size_t >( 0 ) );
    EXPECT_EQ( search.AddFace( IntField{ 1, 2 } ), 0 );
}

TEST( FaceSearchBasicTest, GetFaceNodeRingRejectsOutOfRange )
{
    FaceSearchBasic search;
    EXPECT_THROW( search.GetFaceNodeRing( 0 ), std::out_of_range );
    search.AddFace( IntField{ 1, 2, 3 } );
    EXPECT_NO_THROW( search.GetFaceNodeRing( 0 ) );
    EXPECT_THROW( search.GetFaceNodeRing( 1 ), std::out_of_range );
}

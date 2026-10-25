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

TEST( FaceSearchBasicTest, ClearResetsDerivedSplitMappings )
{
    // Intent: clearing a FaceSearch invalidates mappings derived from its face IDs.
    FaceSearch search;
    ASSERT_EQ( search.AddFace( IntField{ 1, 2, 3, 4 } ), 0 );

    search.interfaceSplit.childFaceIdsOfParent = LinkField{ IntField{ 4, 5 } };
    search.interfaceSplit.childFaceNodeIds = LinkField{ IntField{ 0, 1, 2 } };
    search.interfaceSplit.childFaceNodeFlags = LinkField{ IntField{ 1, -1, 1 } };

    search.Clear();

    EXPECT_EQ( search.Size(), static_cast< std::size_t >( 0 ) );
    EXPECT_TRUE( search.interfaceSplit.childFaceIdsOfParent.empty() );
    EXPECT_TRUE( search.interfaceSplit.childFaceNodeIds.empty() );
    EXPECT_TRUE( search.interfaceSplit.childFaceNodeFlags.empty() );
    EXPECT_EQ( search.AddFace( IntField{ 1, 2 } ), 0 );
}

TEST( FaceSearchBasicTest, RepeatedQuadSplitReplacesParentChildMapping )
{
    // Intent: rebuilding one parent's split mapping must not append duplicate child IDs.
    FaceSearch search;
    ASSERT_EQ( search.AddFace( IntField{ 10, 20, 30, 40 } ), 0 );
    EXPECT_EQ( search.AddFace( IntField{ 10, 20, 30 } ), 1 );
    EXPECT_EQ( search.AddFace( IntField{ 30, 40, 10 } ), 2 );
    EXPECT_EQ( search.AddFace( IntField{ 40, 10, 20 } ), 3 );
    EXPECT_EQ( search.AddFace( IntField{ 20, 30, 40 } ), 4 );
    search.interfaceSplit.ResetForParentFaceCount( search.Size() );

    search.SplitQuad2Tri( 0 );
    search.SplitQuad2Tri( 0 );

    ASSERT_EQ( search.interfaceSplit.childFaceIdsOfParent[ 0 ].size(),
               static_cast< std::size_t >( 4 ) );
    EXPECT_EQ( search.interfaceSplit.childFaceIdsOfParent[ 0 ],
               ( IntField{ 1, 2, 3, 4 } ) );
}

TEST( FaceSearchBasicTest, LineSplitMappingPreservesEdgeDirectionAndEndpointFlags )
{
    // Intent: splitting A->B produces A->M and M->B in both local and global IDs.
    LinkField localLineId{ IntField{ 9, 9 } };
    LinkField localLineFlag{ IntField{ 9, 9 } };
    LinkField lineId{ IntField{ 9, 9 } };

    BuildLineSplitMapping( 12, 18, 7, 25, localLineId, localLineFlag, lineId );

    EXPECT_EQ( localLineId, ( LinkField{ IntField{ 0, 7 }, IntField{ 7, 1 } } ) );
    EXPECT_EQ( localLineFlag, ( LinkField{ IntField{ 1, -1 }, IntField{ -1, 1 } } ) );
    EXPECT_EQ( lineId, ( LinkField{ IntField{ 12, 25 }, IntField{ 25, 18 } } ) );
}

TEST( FaceSearchBasicTest, LineSplitMappingReplacesPreviousOutput )
{
    // Intent: callers may reuse output buffers; stale edges must not survive.
    LinkField localLineId{ IntField{ 99, 98 } };
    LinkField localLineFlag{ IntField{ 97, 96 } };
    LinkField lineId{ IntField{ 95, 94 } };

    BuildLineSplitMapping( 3, 4, 2, 8, localLineId, localLineFlag, lineId );

    ASSERT_EQ( localLineId.size(), static_cast< std::size_t >( 2 ) );
    ASSERT_EQ( localLineFlag.size(), static_cast< std::size_t >( 2 ) );
    ASSERT_EQ( lineId.size(), static_cast< std::size_t >( 2 ) );
    EXPECT_EQ( lineId[ 0 ], ( IntField{ 3, 8 } ) );
    EXPECT_EQ( lineId[ 1 ], ( IntField{ 8, 4 } ) );
}

TEST( FaceSearchBasicTest, GetFaceNodeRingRejectsOutOfRange )
{
    FaceSearchBasic search;
    EXPECT_THROW( search.GetFaceNodeRing( 0 ), std::out_of_range );
    search.AddFace( IntField{ 1, 2, 3 } );
    EXPECT_NO_THROW( search.GetFaceNodeRing( 0 ) );
    EXPECT_THROW( search.GetFaceNodeRing( 1 ), std::out_of_range );
}

TEST( FaceSearchBasicTest, ClearParentMappingOnlyClearsRequestedParent )
{
    InterfaceFaceSplitTable split;
    split.ResetForParentFaceCount( 2 );
    split.childFaceIdsOfParent[ 0 ] = IntField{ 4, 5 };
    split.childFaceIdsOfParent[ 1 ] = IntField{ 7, 8 };

    split.ClearParentMapping( 0 );

    EXPECT_TRUE( split.childFaceIdsOfParent[ 0 ].empty() );
    EXPECT_EQ( split.childFaceIdsOfParent[ 1 ], ( IntField{ 7, 8 } ) );
}

TEST( FaceSearchBasicTest, ClearParentMappingRejectsInvalidParentId )
{
    InterfaceFaceSplitTable split;
    split.ResetForParentFaceCount( 1 );

    EXPECT_THROW( split.ClearParentMapping( -1 ), std::out_of_range );
    EXPECT_THROW( split.ClearParentMapping( 1 ), std::out_of_range );
}

// tests/geometry/QuadDiagonalSplitTest.cpp
// Locks quad->triangle candidate topology used by FaceSearch::SplitQuad2Tri.

#include <gtest/gtest.h>
#include <stdexcept>

#include "FaceSearch.h"
#include "HXDefine.h"

using namespace ONEFLOW;

TEST( QuadDiagonalSplitTest, BuildsFourCandidateTrisForBothDiagonals )
{
    LinkField localTriId, localTriFlag;
    BuildQuadDiagonalLocalTris( localTriId, localTriFlag );

    ASSERT_EQ( localTriId.size(), static_cast< std::size_t >( 4 ) );
    ASSERT_EQ( localTriFlag.size(), static_cast< std::size_t >( 4 ) );

    // Diagonal 0-2
    EXPECT_EQ( localTriId[ 0 ], ( IntField{ 0, 1, 2 } ) );
    EXPECT_EQ( localTriId[ 1 ], ( IntField{ 2, 3, 0 } ) );
    // Diagonal 1-3
    EXPECT_EQ( localTriId[ 2 ], ( IntField{ 3, 0, 1 } ) );
    EXPECT_EQ( localTriId[ 3 ], ( IntField{ 1, 2, 3 } ) );
}

TEST( QuadDiagonalSplitTest, ClearsPreviousContents )
{
    LinkField localTriId, localTriFlag;
    localTriId.push_back( IntField{ 9, 9, 9 } );
    BuildQuadDiagonalLocalTris( localTriId, localTriFlag );
    EXPECT_EQ( localTriId.size(), static_cast< std::size_t >( 4 ) );
    EXPECT_EQ( localTriId[ 0 ], ( IntField{ 0, 1, 2 } ) );
}

TEST( QuadDiagonalSplitTest, MapsLocalCornersThroughQuadNodes )
{
    LinkField localTriId, localTriFlag;
    BuildQuadDiagonalLocalTris( localTriId, localTriFlag );

    // Global nodes on corners 0..3
    const IntField nodeId{ 10, 20, 30, 40 };
    LinkField triId;
    MapLocalTrisToGlobalNodes( nodeId, localTriId, triId );

    ASSERT_EQ( triId.size(), static_cast< std::size_t >( 4 ) );
    EXPECT_EQ( triId[ 0 ], ( IntField{ 10, 20, 30 } ) );
    EXPECT_EQ( triId[ 1 ], ( IntField{ 30, 40, 10 } ) );
    EXPECT_EQ( triId[ 2 ], ( IntField{ 40, 10, 20 } ) );
    EXPECT_EQ( triId[ 3 ], ( IntField{ 20, 30, 40 } ) );
}

TEST( QuadDiagonalSplitTest, MapRejectsBadLocalCorner )
{
    LinkField localTriId;
    localTriId.push_back( IntField{ 0, 1, 5 } ); // 5 out of range for 4-node quad
    const IntField nodeId{ 10, 20, 30, 40 };
    LinkField triId;
    EXPECT_THROW(
        MapLocalTrisToGlobalNodes( nodeId, localTriId, triId ),
        std::out_of_range );
}

TEST( QuadDiagonalSplitTest, MapRejectsShortTriangleWithoutChangingOutput )
{
    LinkField localTriId{ IntField{ 0, 1, 2 }, IntField{ 0, 1 } };
    const IntField nodeId{ 10, 20, 30, 40 };
    LinkField triId{ IntField{ 90, 91, 92 } };

    EXPECT_THROW(
        MapLocalTrisToGlobalNodes( nodeId, localTriId, triId ),
        std::logic_error );

    EXPECT_EQ( triId, ( LinkField{ IntField{ 90, 91, 92 } } ) );
}

TEST( QuadDiagonalSplitTest, ValidatesAllTrianglesBeforeReplacingOutput )
{
    LinkField localTriId{ IntField{ 0, 1, 2 }, IntField{ 1, 2, 9 } };
    const IntField nodeId{ 10, 20, 30, 40 };
    LinkField triId{ IntField{ 90, 91, 92 } };

    EXPECT_THROW(
        MapLocalTrisToGlobalNodes( nodeId, localTriId, triId ),
        std::out_of_range );

    EXPECT_EQ( triId, ( LinkField{ IntField{ 90, 91, 92 } } ) );
}

// tests/geometry/FaceTopoReorderLinkTest.cpp
// ReorderLink must place boundary faces (right cell INVALID) before interior
// faces while preserving relative order within each group.

#include <gtest/gtest.h>

#include "Constant.h"
#include "FaceTopo.h"
#include "HXDefine.h"

using namespace ONEFLOW;

namespace {

void PushFace( FaceTopo & topo, IntField nodes, int left, int right, int type = 5 )
{
    topo.GetFaces().push_back( std::move( nodes ) );
    topo.GetFaceTypes().push_back( type );
    topo.GetLeftCells().push_back( left );
    topo.GetRightCells().push_back( right );
}

} // namespace

TEST( FaceTopoReorderLinkTest, BoundaryFacesMoveToFrontPreservingOrder )
{
    FaceTopo topo;
    // Order: interior, boundary, interior, boundary
    PushFace( topo, IntField{ 0, 1 }, 0, 1 );               // interior A
    PushFace( topo, IntField{ 1, 2 }, 0, INVALID_INDEX );   // boundary B
    PushFace( topo, IntField{ 2, 3 }, 1, 2 );               // interior C
    PushFace( topo, IntField{ 3, 4 }, 2, INVALID_INDEX );   // boundary D

    topo.ReorderLink();

    ASSERT_EQ( topo.GetFaces().size(), static_cast< std::size_t >( 4 ) );
    // Boundary first (B then D), then interior (A then C)
    EXPECT_EQ( topo.GetRightCells()[ 0 ], INVALID_INDEX );
    EXPECT_EQ( topo.GetRightCells()[ 1 ], INVALID_INDEX );
    EXPECT_EQ( topo.GetRightCells()[ 2 ], 1 );
    EXPECT_EQ( topo.GetRightCells()[ 3 ], 2 );

    EXPECT_EQ( topo.GetFaces()[ 0 ], ( IntField{ 1, 2 } ) );
    EXPECT_EQ( topo.GetFaces()[ 1 ], ( IntField{ 3, 4 } ) );
    EXPECT_EQ( topo.GetFaces()[ 2 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( topo.GetFaces()[ 3 ], ( IntField{ 2, 3 } ) );
}

TEST( FaceTopoReorderLinkTest, AllBoundaryUnchangedOrder )
{
    FaceTopo topo;
    PushFace( topo, IntField{ 0, 1 }, 0, INVALID_INDEX );
    PushFace( topo, IntField{ 1, 2 }, 1, INVALID_INDEX );

    topo.ReorderLink();

    EXPECT_EQ( topo.GetFaces()[ 0 ], ( IntField{ 0, 1 } ) );
    EXPECT_EQ( topo.GetFaces()[ 1 ], ( IntField{ 1, 2 } ) );
}

TEST( FaceTopoReorderLinkTest, RejectsInconsistentLeftCellSize )
{
    FaceTopo topo;
    PushFace( topo, IntField{ 0, 1 }, 0, INVALID_INDEX );
    topo.GetLeftCells().pop_back();

    EXPECT_THROW( topo.ReorderLink(), std::logic_error );
}

TEST( FaceTopoReorderLinkTest, ClearsCellToFaceAfterReorder )
{
    FaceTopo topo;
    PushFace( topo, IntField{ 0, 1 }, 0, 1 );
    PushFace( topo, IntField{ 1, 2 }, 0, INVALID_INDEX );
    topo.c2f.push_back( IntField{ 0, 1 } );

    topo.ReorderLink();

    EXPECT_TRUE( topo.c2f.empty() );
}

TEST( BuildBoundaryFirstPermutationTest, OrdersBoundaryThenInterior )
{
    IntField rightCells{ 1, INVALID_INDEX, 2, INVALID_INDEX };
    const IntField perm = BuildBoundaryFirstPermutation( rightCells );

    ASSERT_EQ( perm.size(), static_cast< std::size_t >( 4 ) );
    EXPECT_EQ( perm[ 0 ], 1 );
    EXPECT_EQ( perm[ 1 ], 3 );
    EXPECT_EQ( perm[ 2 ], 0 );
    EXPECT_EQ( perm[ 3 ], 2 );
}

TEST( BuildBoundaryFirstPermutationTest, EmptyInput )
{
    IntField empty;
    const IntField perm = BuildBoundaryFirstPermutation( empty );
    EXPECT_TRUE( perm.empty() );
}

TEST( BuildBoundaryFirstPermutationTest, AllInteriorPreservesOrder )
{
    IntField rightCells{ 1, 2, 3 };
    const IntField perm = BuildBoundaryFirstPermutation( rightCells );
    EXPECT_EQ( perm[ 0 ], 0 );
    EXPECT_EQ( perm[ 1 ], 1 );
    EXPECT_EQ( perm[ 2 ], 2 );
}

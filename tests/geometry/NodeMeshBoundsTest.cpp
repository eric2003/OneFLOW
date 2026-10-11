#include <gtest/gtest.h>
#include "NodeMesh.h"

TEST(NodeMeshBoundsTest, CalculatesAndCachesBounds)
{
    ONEFLOW::NodeMesh mesh;
    mesh.AddPoint( 1.0, -2.0, 3.0 );
    mesh.AddPoint( -4.0, 5.0, 0.0 );

    EXPECT_FALSE( mesh.boxFlag );
    mesh.CalcMinMaxBox();

    EXPECT_TRUE( mesh.boxFlag );
    EXPECT_DOUBLE_EQ( mesh.pmin[ 0 ], -4.0 );
    EXPECT_DOUBLE_EQ( mesh.pmin[ 1 ], -2.0 );
    EXPECT_DOUBLE_EQ( mesh.pmin[ 2 ], 0.0 );
    EXPECT_DOUBLE_EQ( mesh.pmax[ 0 ], 1.0 );
    EXPECT_DOUBLE_EQ( mesh.pmax[ 1 ], 5.0 );
    EXPECT_DOUBLE_EQ( mesh.pmax[ 2 ], 3.0 );
}

TEST(NodeMeshBoundsTest, AddingPointInvalidatesCachedBounds)
{
    ONEFLOW::NodeMesh mesh;
    mesh.AddPoint( 0.0, 0.0, 0.0 );
    mesh.CalcMinMaxBox();
    ASSERT_TRUE( mesh.boxFlag );

    mesh.AddPoint( 2.0, 3.0, 4.0 );
    EXPECT_FALSE( mesh.boxFlag );

    mesh.CalcMinMaxBox();
    EXPECT_DOUBLE_EQ( mesh.pmax[ 0 ], 2.0 );
    EXPECT_DOUBLE_EQ( mesh.pmax[ 1 ], 3.0 );
    EXPECT_DOUBLE_EQ( mesh.pmax[ 2 ], 4.0 );
}

TEST(NodeMeshBoundsTest, EmptyMeshHasDeterministicBounds)
{
    ONEFLOW::NodeMesh mesh;
    mesh.CalcMinMaxBox();

    ASSERT_TRUE( mesh.boxFlag );
    for ( int iDim = 0; iDim < 3; ++ iDim )
    {
        EXPECT_DOUBLE_EQ( mesh.pmin[ iDim ], 0.0 );
        EXPECT_DOUBLE_EQ( mesh.pmax[ iDim ], 0.0 );
    }
}

TEST(NodeMeshBoundsTest, ResizingInvalidatesCachedBounds)
{
    ONEFLOW::NodeMesh mesh;
    mesh.AddPoint( 1.0, 2.0, 3.0 );
    mesh.CalcMinMaxBox();
    ASSERT_TRUE( mesh.boxFlag );

    mesh.CreateNodes( 2 );
    EXPECT_FALSE( mesh.boxFlag );
}

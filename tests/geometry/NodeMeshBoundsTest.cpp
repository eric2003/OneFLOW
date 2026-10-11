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

TEST(NodeMeshBoundsTest, CopyOwnsIndependentCoordinateStorage)
{
    ONEFLOW::NodeMesh original;
    original.AddPoint( 1.0, 2.0, 3.0 );
    original.CalcMinMaxBox();

    ONEFLOW::NodeMesh copy = original;
    copy.SetPoint( 0, -4.0, -5.0, -6.0 );

    EXPECT_DOUBLE_EQ( original.xN[ 0 ], 1.0 );
    EXPECT_DOUBLE_EQ( original.yN[ 0 ], 2.0 );
    EXPECT_DOUBLE_EQ( original.zN[ 0 ], 3.0 );
    EXPECT_DOUBLE_EQ( copy.xN[ 0 ], -4.0 );
    EXPECT_FALSE( copy.boxFlag );
    EXPECT_TRUE( original.boxFlag );
}

TEST(NodeMeshBoundsTest, PointValueInterfaceReadsSoAStorage)
{
    ONEFLOW::NodeMesh mesh;
    mesh.AddPoint( 1.25, -2.5, 3.75 );
    mesh.AddPoint( 4.0, 5.0, 6.0 );

    const ONEFLOW::NodeMesh & constMesh = mesh;
    const ONEFLOW::NodeMesh::Point3 point = constMesh.GetPoint( 1 );

    EXPECT_DOUBLE_EQ( point.x, 4.0 );
    EXPECT_DOUBLE_EQ( point.y, 5.0 );
    EXPECT_DOUBLE_EQ( point.z, 6.0 );
}

TEST(NodeMeshBoundsTest, SettingPointUpdatesSoAStorageAndInvalidatesBounds)
{
    ONEFLOW::NodeMesh mesh;
    mesh.AddPoint( 1.0, 2.0, 3.0 );
    mesh.CalcMinMaxBox();
    ASSERT_TRUE( mesh.boxFlag );

    mesh.SetPoint( 0, -4.0, 5.0, -6.0 );

    EXPECT_FALSE( mesh.boxFlag );
    EXPECT_DOUBLE_EQ( mesh.xN[ 0 ], -4.0 );
    EXPECT_DOUBLE_EQ( mesh.yN[ 0 ], 5.0 );
    EXPECT_DOUBLE_EQ( mesh.zN[ 0 ], -6.0 );

    mesh.CalcMinMaxBox();
    EXPECT_DOUBLE_EQ( mesh.pmin[ 0 ], -4.0 );
    EXPECT_DOUBLE_EQ( mesh.pmax[ 1 ], 5.0 );
    EXPECT_DOUBLE_EQ( mesh.pmin[ 2 ], -6.0 );
}

TEST(NodeMeshBoundsTest, PointValueInterfaceRejectsOutOfRangeIndex)
{
    ONEFLOW::NodeMesh mesh;
    mesh.AddPoint( 1.0, 2.0, 3.0 );

    EXPECT_THROW( mesh.GetPoint( 1 ), std::out_of_range );
    EXPECT_THROW( mesh.SetPoint( 1, 4.0, 5.0, 6.0 ), std::out_of_range );
}


TEST(NodeMeshBoundsTest, ScaleAndTranslateUpdatesCoordinatesAndInvalidatesBounds)
{
    ONEFLOW::NodeMesh mesh;
    mesh.AddPoint( 1.0, -2.0, 3.0 );
    mesh.CalcMinMaxBox();
    ASSERT_TRUE( mesh.boxFlag );

    mesh.ScaleAndTranslate( 2.0, 10.0, 20.0, -4.0 );

    EXPECT_FALSE( mesh.boxFlag );
    const auto point = mesh.GetPoint( 0 );
    EXPECT_DOUBLE_EQ( point.x, 12.0 );
    EXPECT_DOUBLE_EQ( point.y, 16.0 );
    EXPECT_DOUBLE_EQ( point.z, 2.0 );

    mesh.CalcMinMaxBox();
    EXPECT_DOUBLE_EQ( mesh.pmin[ 0 ], 12.0 );
    EXPECT_DOUBLE_EQ( mesh.pmax[ 1 ], 16.0 );
    EXPECT_DOUBLE_EQ( mesh.pmin[ 2 ], 2.0 );
}

TEST(NodeMeshBoundsTest, AxisTransformUpdatesCoordinatesAndInvalidatesBounds)
{
    ONEFLOW::NodeMesh mesh;
    mesh.AddPoint( 1.0, 2.0, 3.0 );
    mesh.CalcMinMaxBox();
    ASSERT_TRUE( mesh.boxFlag );

    mesh.TurnZAxisToYAxis();

    EXPECT_FALSE( mesh.boxFlag );
    const auto point = mesh.GetPoint( 0 );
    EXPECT_DOUBLE_EQ( point.x, 1.0 );
    EXPECT_DOUBLE_EQ( point.y, 3.0 );
    EXPECT_DOUBLE_EQ( point.z, -2.0 );
}

TEST(NodeMeshBoundsTest, BoundsRejectInconsistentCoordinateArraysEvenWhenCached)
{
    ONEFLOW::NodeMesh mesh;
    mesh.AddPoint( 1.0, 2.0, 3.0 );
    mesh.CalcMinMaxBox();
    ASSERT_TRUE( mesh.boxFlag );

    mesh.yN.clear();

    EXPECT_THROW( mesh.CalcMinMaxBox(), std::logic_error );
}

TEST(NodeMeshBoundsTest, CoordinateTransformsRejectInconsistentArraysBeforeMutation)
{
    ONEFLOW::NodeMesh mesh;
    mesh.AddPoint( 1.0, 2.0, 3.0 );
    mesh.yN.clear();

    EXPECT_THROW(
        mesh.ScaleAndTranslate( 2.0, 10.0, 20.0, 30.0 ),
        std::logic_error );
    EXPECT_DOUBLE_EQ( mesh.xN[ 0 ], 1.0 );

    EXPECT_THROW( mesh.TurnZAxisToYAxis(), std::logic_error );
    EXPECT_DOUBLE_EQ( mesh.xN[ 0 ], 1.0 );
}

TEST(NodeMeshBoundsTest, CreateNodesRejectsNegativeCountWithoutChangingMesh)
{
    ONEFLOW::NodeMesh mesh;
    mesh.AddPoint( 1.0, 2.0, 3.0 );

    EXPECT_THROW( mesh.CreateNodes( -1 ), std::invalid_argument );

    ASSERT_EQ( mesh.GetNumberOfNodes(), 1 );
    EXPECT_DOUBLE_EQ( mesh.xN[ 0 ], 1.0 );
    EXPECT_DOUBLE_EQ( mesh.yN[ 0 ], 2.0 );
    EXPECT_DOUBLE_EQ( mesh.zN[ 0 ], 3.0 );
}

TEST(NodeMeshBoundsTest, SetCoordinatesReplacesStorageAndInvalidatesBounds)
{
    ONEFLOW::NodeMesh mesh;
    mesh.AddPoint( 1.0, 2.0, 3.0 );
    mesh.CalcMinMaxBox();
    ASSERT_TRUE( mesh.boxFlag );

    ONEFLOW::RealField x = { -1.0, -2.0 };
    ONEFLOW::RealField y = { 4.0, 5.0 };
    ONEFLOW::RealField z = { 6.0, 7.0 };
    mesh.SetCoordinates( x, y, z );

    EXPECT_FALSE( mesh.boxFlag );
    ASSERT_EQ( mesh.GetNumberOfNodes(), 2 );
    EXPECT_DOUBLE_EQ( mesh.xN[ 0 ], -1.0 );
    EXPECT_DOUBLE_EQ( mesh.yN[ 1 ], 5.0 );
    EXPECT_DOUBLE_EQ( mesh.zN[ 1 ], 7.0 );

    mesh.CalcMinMaxBox();
    EXPECT_DOUBLE_EQ( mesh.pmin[ 0 ], -2.0 );
    EXPECT_DOUBLE_EQ( mesh.pmax[ 2 ], 7.0 );
}

TEST(NodeMeshBoundsTest, SetCoordinatesRejectsMismatchedSizesWithoutMutation)
{
    ONEFLOW::NodeMesh mesh;
    mesh.AddPoint( 1.0, 2.0, 3.0 );

    ONEFLOW::RealField x = { 4.0, 5.0 };
    ONEFLOW::RealField y = { 6.0 };
    ONEFLOW::RealField z = { 7.0, 8.0 };

    EXPECT_THROW( mesh.SetCoordinates( x, y, z ), std::invalid_argument );
    ASSERT_EQ( mesh.GetNumberOfNodes(), 1 );
    EXPECT_DOUBLE_EQ( mesh.xN[ 0 ], 1.0 );
    EXPECT_DOUBLE_EQ( mesh.yN[ 0 ], 2.0 );
    EXPECT_DOUBLE_EQ( mesh.zN[ 0 ], 3.0 );
}

TEST(NodeMeshBoundsTest, RecalculatesBoundsAfterDirectCoordinateMutation)
{
    ONEFLOW::NodeMesh mesh;
    mesh.AddPoint( 1.0, 2.0, 3.0 );
    mesh.CalcMinMaxBox();
    ASSERT_TRUE( mesh.boxFlag );

    // Public SoA arrays remain part of the existing API, so callers can bypass
    // SetPoint and its cache invalidation.
    mesh.xN[ 0 ] = -8.0;

    mesh.CalcMinMaxBox();

    EXPECT_DOUBLE_EQ( mesh.pmin[ 0 ], -8.0 );
    EXPECT_DOUBLE_EQ( mesh.pmax[ 0 ], -8.0 );
}

TEST(NodeMeshBoundsTest, AddPointRejectsInconsistentCoordinateArraysWithoutMutation)
{
    ONEFLOW::NodeMesh mesh;
    mesh.AddPoint( 1.0, 2.0, 3.0 );
    mesh.yN.clear();

    EXPECT_THROW( mesh.AddPoint( 4.0, 5.0, 6.0 ), std::logic_error );

    ASSERT_EQ( mesh.xN.size(), 1 );
    EXPECT_TRUE( mesh.yN.empty() );
    ASSERT_EQ( mesh.zN.size(), 1 );
    EXPECT_DOUBLE_EQ( mesh.xN[ 0 ], 1.0 );
    EXPECT_DOUBLE_EQ( mesh.zN[ 0 ], 3.0 );
}

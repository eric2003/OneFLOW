#include "GridLayout.h"
#include "GridLayoutParser.h"

#include <gtest/gtest.h>

#include <string>
#include <stdexcept>

namespace ONEFLOW
{

TEST( GridLayoutParserTest, ParsesLaminarPlateLayout )
{
    const std::string fileName =
        std::string( ONEFLOW_SOURCE_DIR ) +
        "/tests/geometry/data/laminarPlate2dLayout2.txt";

    const GridLayout layout = GridLayoutParser().Parse( fileName );

    EXPECT_EQ( layout.points.size(), 6u );
    EXPECT_EQ( layout.lines.size(), 7u );
    EXPECT_EQ( layout.circles.size(), 0u );
    EXPECT_EQ( layout.dimensions.size(), 7u );
    EXPECT_EQ( layout.distributions.size(), 7u );
    EXPECT_EQ( layout.boundaries.size(), 7u );
    EXPECT_EQ( layout.lineToFaces.size(), 8u );
    EXPECT_EQ( layout.faceToBlocks.size(), 2u );

    EXPECT_EQ( layout.points[ 0 ].id, 1 );
    EXPECT_DOUBLE_EQ( layout.points[ 0 ].x, -2.0 );
    EXPECT_EQ( layout.dimensions[ 2 ].pointCount, 85 );
    EXPECT_EQ( layout.distributions[ 0 ].type, GridDistributionType::Distance );
    EXPECT_DOUBLE_EQ( layout.distributions[ 0 ].endValue, 0.01 );
    EXPECT_EQ( layout.boundaries[ 0 ].boundaryType, 3 );
    EXPECT_EQ(
        GridLayout::ClassifyBoundaryType( layout.boundaries[ 0 ].boundaryType ),
        GridBoundaryKind::Physical );
    EXPECT_EQ(
        GridLayout::ClassifyBoundaryType( layout.boundaries[ 3 ].boundaryType ),
        GridBoundaryKind::Interface );
    EXPECT_EQ( layout.lineToFaces[ 0 ].lineId, 1 );
    EXPECT_EQ( layout.faceToBlocks[ 1 ].blockId, 2 );
}

TEST( GridLayoutParserTest, ClassifiesBoundaryTypes )
{
    EXPECT_EQ(
        GridLayout::ClassifyBoundaryType( 0 ),
        GridBoundaryKind::None );
    EXPECT_EQ(
        GridLayout::ClassifyBoundaryType( -1 ),
        GridBoundaryKind::Interface );
    EXPECT_EQ(
        GridLayout::ClassifyBoundaryType( 100 ),
        GridBoundaryKind::Interior );
    EXPECT_EQ(
        GridLayout::ClassifyBoundaryType( 1000 ),
        GridBoundaryKind::Overset );
    EXPECT_EQ(
        GridLayout::ClassifyBoundaryType( 2 ),
        GridBoundaryKind::Physical );
    EXPECT_EQ(
        GridLayout::ClassifyBoundaryType( 61 ),
        GridBoundaryKind::Physical );
    EXPECT_EQ(
        GridLayout::ClassifyBoundaryType( -2 ),
        GridBoundaryKind::Physical );
}

TEST( GridLayoutParserTest, ParsesTwoConnectedThreeDimensionalCubes )
{
    const std::string fileName =
        std::string( ONEFLOW_SOURCE_DIR ) +
        "/grid/two_connected_cubes/grid/twoConnectedCubesLayout.txt";

    const GridLayout layout = GridLayoutParser().Parse( fileName );

    EXPECT_EQ( layout.points.size(), 12u );
    EXPECT_EQ( layout.lines.size(), 20u );
    EXPECT_EQ( layout.dimensions.size(), 20u );
    EXPECT_EQ( layout.distributions.size(), 20u );
    EXPECT_EQ( layout.boundaries.size(), 20u );
    EXPECT_EQ( layout.lineToFaces.size(), 44u );
    EXPECT_EQ( layout.faceToBlocks.size(), 12u );

    int sharedFaceRelations = 0;
    bool firstBlockUsesSharedFace = false;
    bool secondBlockUsesSharedFace = false;
    for ( const auto & relation : layout.faceToBlocks )
    {
        if ( relation.faceId != 4 ) continue;

        ++ sharedFaceRelations;
        firstBlockUsesSharedFace =
            firstBlockUsesSharedFace ||
            ( relation.blockId == 1 && relation.position == 2 );
        secondBlockUsesSharedFace =
            secondBlockUsesSharedFace ||
            ( relation.blockId == 2 && relation.position == 1 );
    }
    EXPECT_EQ( sharedFaceRelations, 2 );
    EXPECT_TRUE( firstBlockUsesSharedFace );
    EXPECT_TRUE( secondBlockUsesSharedFace );
}

TEST( GridLayoutParserTest, ParsesLaminarFlatPlateLayout )
{
    const std::string fileName =
        std::string( ONEFLOW_SOURCE_DIR ) +
        "/grid/laminar_flat_plate/grid/laminarFlatPlateLayout.txt";

    const GridLayout layout = GridLayoutParser().Parse( fileName );

    EXPECT_EQ( layout.points.size(), 6u );
    EXPECT_EQ( layout.lines.size(), 7u );
    EXPECT_EQ( layout.dimensions.size(), 7u );
    EXPECT_EQ( layout.distributions.size(), 7u );
    EXPECT_EQ( layout.boundaries.size(), 7u );
    EXPECT_EQ( layout.lineToFaces.size(), 8u );
    EXPECT_EQ( layout.faceToBlocks.size(), 2u );

    EXPECT_EQ( layout.dimensions[ 0 ].pointCount, 33 );
    EXPECT_EQ( layout.dimensions[ 1 ].pointCount, 41 );
    EXPECT_EQ( layout.dimensions[ 2 ].pointCount, 85 );
    EXPECT_EQ( layout.boundaries[ 0 ].boundaryType, 3 );
    EXPECT_EQ( layout.boundaries[ 1 ].boundaryType, 2 );
    EXPECT_EQ( layout.boundaries[ 3 ].boundaryType, -1 );
    EXPECT_EQ( layout.boundaries[ 4 ].boundaryType, 6 );

    EXPECT_EQ(
        GridLayout::ClassifyBoundaryType( layout.boundaries[ 1 ].boundaryType ),
        GridBoundaryKind::Physical );
    EXPECT_EQ(
        GridLayout::ClassifyBoundaryType( layout.boundaries[ 3 ].boundaryType ),
        GridBoundaryKind::Interface );
}

TEST( GridLayoutParserTest, RejectsUnknownPointReference )
{
    const std::string fileName =
        std::string( ONEFLOW_SOURCE_DIR ) +
        "/tests/geometry/data/laminarPlate2dLayoutInvalid.txt";

    EXPECT_THROW(
        GridLayoutParser().Parse( fileName ),
        std::runtime_error );
}

} // namespace ONEFLOW

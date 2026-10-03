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
    EXPECT_DOUBLE_EQ( layout.distributions[ 0 ].value2, 0.01 );
    EXPECT_EQ( layout.boundaries[ 0 ].boundaryType, 3 );
    EXPECT_EQ( layout.lineToFaces[ 0 ].lineId, 1 );
    EXPECT_EQ( layout.faceToBlocks[ 1 ].blockId, 2 );
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

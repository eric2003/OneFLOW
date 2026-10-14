#include <gtest/gtest.h>

#include <sstream>
#include <string>

#include "GridSummary.h"
#include "GridSummaryJson.h"

TEST( GridSummaryTest, AcceptsValidOneDimensionalSummary )
{
    const ONEFLOW::GridSummary summary{ 1, 4, 5, 3, 2 };

    EXPECT_TRUE( ONEFLOW::IsValidGridSummary( summary ) );
    EXPECT_TRUE( ONEFLOW::ValidateGridSummary( summary ).empty() );
}

TEST( GridSummaryTest, AcceptsValidTwoAndThreeDimensionalSummaries )
{
    const ONEFLOW::GridSummary twoDimensional{ 2, 12, 17, 6, 12 };
    const ONEFLOW::GridSummary threeDimensional{ 3, 8, 12, 5, 6 };

    EXPECT_TRUE( ONEFLOW::IsValidGridSummary( twoDimensional ) );
    EXPECT_TRUE( ONEFLOW::IsValidGridSummary( threeDimensional ) );
}

TEST( GridSummaryTest, ReportsEveryInvalidCount )
{
    const ONEFLOW::GridSummary summary{ 4, -1, -2, -3, -4 };

    const auto errors = ONEFLOW::ValidateGridSummary( summary );

    EXPECT_EQ( errors.size(), 5u );
    EXPECT_FALSE( ONEFLOW::IsValidGridSummary( summary ) );
}

TEST( GridSummaryTest, RejectsBoundaryFaceCountGreaterThanFaceCount )
{
    const ONEFLOW::GridSummary summary{ 2, 4, 3, 1, 4 };

    const auto errors = ONEFLOW::ValidateGridSummary( summary );

    ASSERT_EQ( errors.size(), 1u );
    EXPECT_EQ( errors[ 0 ], "boundary face count must not exceed face count" );
}

TEST( GridSummaryTest, RejectsFacesOrCellsWithoutNodes )
{
    const ONEFLOW::GridSummary summary{ 1, 0, 1, 1, 1 };

    const auto errors = ONEFLOW::ValidateGridSummary( summary );

    ASSERT_EQ( errors.size(), 1u );
    EXPECT_EQ( errors[ 0 ], "faces or cells require at least one node" );
}

TEST( GridSummaryTest, AllowsAnEmptyGridWithValidDimension )
{
    const ONEFLOW::GridSummary summary{ 3, 0, 0, 0, 0 };

    EXPECT_TRUE( ONEFLOW::IsValidGridSummary( summary ) );
}

TEST( GridSummaryJsonTest, WritesVersionedMachineReadableMetadata )
{
    const ONEFLOW::GridSummary summary{ 2, 12, 17, 6, 12 };
    std::ostringstream output;

    ONEFLOW::WriteGridSummaryJson( output, summary );

    EXPECT_EQ(
        output.str(),
        "{\n"
        "  \"schema\": \"oneflow.grid-summary\",\n"
        "  \"version\": 1,\n"
        "  \"dimension\": 2,\n"
        "  \"counts\": {\n"
        "    \"nodes\": 12,\n"
        "    \"faces\": 17,\n"
        "    \"cells\": 6,\n"
        "    \"boundaryFaces\": 12\n"
        "  }\n"
        "}\n" );
}

TEST( GridSummaryJsonTest, SupportsWritingToAnExistingStream )
{
    const ONEFLOW::GridSummary summary{ 1, 4, 5, 3, 2 };
    std::ostringstream output;
    output << "prefix:";

    ONEFLOW::WriteGridSummaryJson( output, summary );

    EXPECT_EQ( output.str().find( "prefix:" ), 0u );
    EXPECT_NE( output.str().find( "\"schema\": \"oneflow.grid-summary\"" ),
               std::string::npos );
}

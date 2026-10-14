#include <gtest/gtest.h>

#include <sstream>
#include <stdexcept>
#include <string>

#include <nlohmann/json.hpp>

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

    const nlohmann::json document = nlohmann::json::parse( output.str() );
    EXPECT_EQ( document.at( "schema" ).get< std::string >(), "oneflow.grid-summary" );
    EXPECT_EQ( document.at( "version" ).get< int >(), 1 );
    EXPECT_EQ( document.at( "dimension" ).get< int >(), 2 );
    EXPECT_EQ( document.at( "counts" ).at( "nodes" ).get< int >(), 12 );
    EXPECT_EQ( document.at( "counts" ).at( "faces" ).get< int >(), 17 );
    EXPECT_EQ( document.at( "counts" ).at( "cells" ).get< int >(), 6 );
    EXPECT_EQ( document.at( "counts" ).at( "boundaryFaces" ).get< int >(), 12 );
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

TEST( GridSummaryJsonTest, RejectsInvalidSummaryBeforeWriting )
{
    const ONEFLOW::GridSummary summary{ 4, 1, 1, 1, 0 };
    std::ostringstream output;
    output << "prefix:";

    EXPECT_THROW(
        ONEFLOW::WriteGridSummaryJson( output, summary ),
        std::invalid_argument );
    EXPECT_EQ( output.str(), "prefix:" );
}

#include <gtest/gtest.h>
#include "FieldSimu.h"

using ONEFLOW::FieldPipeline;

TEST( FieldPipelineMeta, StageCountIsSix )
{
    EXPECT_EQ( FieldPipeline::kStageCount, 6 );
}

TEST( FieldPipelineMeta, StageNamesMatchPipelineOrder )
{
    EXPECT_STREQ( FieldPipeline::StageName( 0 ), "SetupGlobals" );
    EXPECT_STREQ( FieldPipeline::StageName( 1 ), "LoadGrid" );
    EXPECT_STREQ( FieldPipeline::StageName( 2 ), "PrepareWallDist" );
    EXPECT_STREQ( FieldPipeline::StageName( 3 ), "CreateSolvers" );
    EXPECT_STREQ( FieldPipeline::StageName( 4 ), "InitFlowField" );
    EXPECT_STREQ( FieldPipeline::StageName( 5 ), "Run" );
}

TEST( FieldPipelineMeta, OutOfRangeStageNameIsEmpty )
{
    EXPECT_STREQ( FieldPipeline::StageName( -1 ), "" );
    EXPECT_STREQ( FieldPipeline::StageName( 6 ), "" );
}
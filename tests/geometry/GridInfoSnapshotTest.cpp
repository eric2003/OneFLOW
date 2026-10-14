#include <gtest/gtest.h>

#include "GridInfoJsonWriter.h"
#include "GridInfoSnapshotAdapter.h"

#include <nlohmann/json.hpp>

#include <sstream>

TEST( GridInfoSnapshotTest, EmptyRuntimeCollectionProducesEmptySnapshot )
{
    const ONEFLOW::Grids grids;
    const ONEFLOW::GridInfoSnapshot snapshot = ONEFLOW::MakeGridInfoSnapshot( grids );

    EXPECT_EQ( snapshot.dimension, 0 );
    EXPECT_EQ( snapshot.nodeCount, 0 );
    EXPECT_EQ( snapshot.faceCount, 0 );
    EXPECT_EQ( snapshot.cellCount, 0 );
    EXPECT_EQ( snapshot.boundaryFaceCount, 0 );
    EXPECT_FALSE( snapshot.hasBoundingBox );
    EXPECT_TRUE( snapshot.zones.empty() );
}

TEST( GridInfoSnapshotTest, JsonWriterSerializesSnapshotWithoutRuntimeGridAccess )
{
    ONEFLOW::GridInfoSnapshot snapshot;
    snapshot.dimension = 2;
    snapshot.nodeCount = 4;
    snapshot.faceCount = 4;
    snapshot.cellCount = 1;
    snapshot.boundaryFaceCount = 4;
    snapshot.hasBoundingBox = true;
    snapshot.minimum = { 0.0, -1.0, 0.0 };
    snapshot.maximum = { 2.0, 1.0, 0.0 };
    snapshot.zones.push_back( { "plate", 2, 4, 4, 1, 4 } );

    std::ostringstream output;
    ONEFLOW::WriteGridInfoJson( output, snapshot );
    const nlohmann::json document = nlohmann::json::parse( output.str() );

    EXPECT_EQ( document.at( "schema" ), "oneflow.mesh-info" );
    EXPECT_EQ( document.at( "version" ), 1 );
    EXPECT_EQ( document.at( "dimension" ), 2 );
    EXPECT_EQ( document.at( "counts" ).at( "zones" ), 1 );
    EXPECT_EQ( document.at( "counts" ).at( "nodes" ), 4 );
    EXPECT_EQ( document.at( "boundingBox" ).at( "min" ),
               nlohmann::json( { 0.0, -1.0, 0.0 } ) );
    EXPECT_EQ( document.at( "zones" ).at( 0 ).at( "name" ), "plate" );
}

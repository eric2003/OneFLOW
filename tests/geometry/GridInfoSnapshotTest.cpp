#include <gtest/gtest.h>

#include "GridInfoJsonWriter.h"
#include "GridInfoSnapshotAdapter.h"

#include <nlohmann/json.hpp>

#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include "UnsGrid.h"
#include "NodeMesh.h"

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
    snapshot.caseRoot = "/project/cases/mesh";
    snapshot.projectRoot = "/project";
    snapshot.inputRoot = "/project/cases/source";
    snapshot.outputRoot = "/project/results";
    snapshot.sourceGridPath = "/project/cases/source/grid.grd";
    snapshot.outputGridPath = "/project/results/grid.ofl";
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
    EXPECT_EQ( document.at( "paths" ).at( "inputRoot" ), "/project/cases/source" );
    EXPECT_EQ( document.at( "paths" ).at( "outputRoot" ), "/project/results" );
    EXPECT_EQ( document.at( "paths" ).at( "sourceGrid" ), "/project/cases/source/grid.grd" );
    EXPECT_EQ( document.at( "paths" ).at( "outputGrid" ), "/project/results/grid.ofl" );
    EXPECT_EQ( document.at( "counts" ).at( "zones" ), 1 );
    EXPECT_EQ( document.at( "counts" ).at( "nodes" ), 4 );
    EXPECT_EQ( document.at( "boundingBox" ).at( "min" ),
               nlohmann::json( { 0.0, -1.0, 0.0 } ) );
    EXPECT_EQ( document.at( "zones" ).at( 0 ).at( "name" ), "plate" );
}

TEST( GridInfoSnapshotTest, RejectsNullZoneInsteadOfDereferencingIt )
{
    ONEFLOW::Grids grids;
    grids.push_back( nullptr );

    EXPECT_THROW( ONEFLOW::MakeGridInfoSnapshot( grids ), std::invalid_argument );
}

TEST( GridInfoSnapshotTest, RejectsNonFiniteNodeCoordinates )
{
    ONEFLOW::Grids grids;
    auto grid = std::make_unique< ONEFLOW::UnsGrid >();
    grid->Init();
    grid->dimension = 2;
    grid->nodeMesh->CreateNodes( 1 );
    grid->nodeMesh->SetPoint(
        0, std::numeric_limits< ONEFLOW::Real >::quiet_NaN(), 0.0, 0.0 );
    grids.push_back( std::move( grid ) );

    EXPECT_THROW( ONEFLOW::MakeGridInfoSnapshot( grids ), std::domain_error );
}

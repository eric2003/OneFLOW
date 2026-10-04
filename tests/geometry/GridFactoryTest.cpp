// tests/geometry/GridFactoryTest.cpp
// Lifecycle + typed config / mediator unit tests (no file I/O, no global DataBase).

#include <gtest/gtest.h>

#include "GridFactory.h"
#include "GridMediator.h"
#include "GridTypes.h"

#include <memory>
#include <stdexcept>
#include <string>
#include <set>

using namespace ONEFLOW;

// ---------------------------------------------------------------------------
// GridFactory lifecycle (existing characterization)
// ---------------------------------------------------------------------------

TEST( GridFactoryLifecycleTest, StackAllocationIsSafe )
{
    EXPECT_NO_THROW( {
        GridFactory gf;
        // Do not call Run(): avoids DataBase / file I/O in unit tests.
    } );
}

TEST( GridFactoryLifecycleTest, HeapAllocationAlsoWorksButIsDiscouraged )
{
    EXPECT_NO_THROW( {
        auto * gf = new GridFactory();
        delete gf;
    } );
}

// ---------------------------------------------------------------------------
// GridTypes: objective / file type parsing (header-only helpers)
// ---------------------------------------------------------------------------

TEST( GridTypesTest, ParseGridObjectiveKnownValues )
{
    EXPECT_EQ( ParseGridObjective( 0 ), GridObjective::GenerateClassic );
    EXPECT_EQ( ParseGridObjective( 1 ), GridObjective::ConvertOnly );
    EXPECT_EQ( ParseGridObjective( 2 ), GridObjective::GenerateInp );
    EXPECT_EQ( ParseGridObjective( 3 ), GridObjective::Partition );
}

TEST( GridTypesTest, ParseGridObjectiveUnknownIsNullopt )
{
    EXPECT_FALSE( ParseGridObjective( -1 ).has_value() );
    EXPECT_FALSE( ParseGridObjective( 99 ).has_value() );
}

TEST( GridTypesTest, ParseGridFileTypeCaseInsensitive )
{
    EXPECT_EQ( ParseGridFileType( "plot3d" ), GridFileType::Plot3D );
    EXPECT_EQ( ParseGridFileType( "PLOT3D" ), GridFileType::Plot3D );
    EXPECT_EQ( ParseGridFileType( "su2" ), GridFileType::SU2 );
    EXPECT_EQ( ParseGridFileType( "CGNS" ), GridFileType::CGNS );
    EXPECT_EQ( ParseGridFileType( "oneflow" ), GridFileType::OneFLOW );
    EXPECT_EQ( ParseGridFileType( "not-a-format" ), GridFileType::Unknown );
}

TEST( GridTypesTest, ToStringRoundTripTokens )
{
    EXPECT_EQ( ToString( GridFileType::Plot3D ), "plot3d" );
    EXPECT_EQ( ToString( GridFileType::SU2 ), "su2" );
    EXPECT_EQ( ToString( GridFileType::CGNS ), "cgns" );
    EXPECT_EQ( ToString( GridFileType::OneFLOW ), "oneflow" );

    EXPECT_EQ( ToString( GridObjective::GenerateClassic ), "GenerateClassic" );
    EXPECT_EQ( ToString( GridOp::CalcMetrics ), "CALC_METRICS" );
}

TEST( GridTypesTest, ParseGridOpKnownTokens )
{
    EXPECT_EQ( ParseGridOp( "CALC_METRICS" ), GridOp::CalcMetrics );
    EXPECT_EQ( ParseGridOp( "fill_wall_struct" ), GridOp::FillWallStruct );
    EXPECT_FALSE( ParseGridOp( "NOT_A_REAL_OP" ).has_value() );
}

TEST( GridTypesTest, GridConfigDefaultValues )
{
    GridConfig cfg;
    EXPECT_EQ( cfg.objective, GridObjective::ConvertOnly );
    EXPECT_EQ( cfg.sourceType, GridFileType::Unknown );
    EXPECT_EQ( cfg.targetType, GridFileType::Unknown );
    EXPECT_EQ( cfg.assemblyMode, GridAssemblyMode::AggregateZones );
    EXPECT_EQ( cfg.topology, GridTopology::Unknown );
    EXPECT_EQ( cfg.axisDirection, GridAxisDirection::Y );
    EXPECT_EQ( cfg.scale, 1.0 );
    EXPECT_TRUE( cfg.sourceFile.empty() );
}

// ---------------------------------------------------------------------------
// ZgridMediator RAII (no file I/O)
// ---------------------------------------------------------------------------

TEST( ZgridMediatorTest, CreateSimpleOwnsMediator )
{
    ZgridMediator zgm;
    EXPECT_EQ( zgm.GetSize(), 0 );

    zgm.CreateSimple( 3 );
    ASSERT_EQ( zgm.GetSize(), 1 );

    GridMediator * gm = &zgm.GetGridMediator( 0 );
    ASSERT_NE( gm, nullptr );
    EXPECT_EQ( gm->numberOfZones, 3 );
    // Destructor of zgm must free the unique_ptr without leak/crash.
}

TEST( ZgridMediatorTest, AddUniquePtrTakesOwnership )
{
    ZgridMediator zgm;
    auto owned = std::make_unique< GridMediator >();
    owned->numberOfZones = 7;
    owned->gridType = "plot3d";

    zgm.AddGridMediator( std::move( owned ) );
    EXPECT_EQ( owned, nullptr );
    ASSERT_EQ( zgm.GetSize(), 1 );
    EXPECT_EQ( zgm.GetGridMediator( 0 )->numberOfZones, 7 );
    EXPECT_EQ( zgm.GetGridMediator( 0 )->gridType, "plot3d" );
}

TEST( ZgridMediatorTest, AddUniquePtrTakesOwnershipWithTwoZones )
{
    ZgridMediator zgm;
    auto owned = std::make_unique< GridMediator >();
    owned->numberOfZones = 2;
    zgm.AddGridMediator( std::move( owned ) );
    EXPECT_EQ( owned, nullptr );
    ASSERT_EQ( zgm.GetSize(), 1 );
    EXPECT_EQ( zgm.GetGridMediator( 0 ).numberOfZones, 2 );
}

TEST( GridFactoryDispatchTest, UnknownObjectiveThrows )
{
    GridFactory gf;
    GridConfig cfg;
    // Force an out-of-range objective without going through ParseGridObjective.
    cfg.objective = static_cast< GridObjective >( 42 );

    EXPECT_THROW( gf.Run( cfg ), std::invalid_argument );
}


TEST( GridOpCatalogTest, AllTokensRoundTrip )
{
    for ( std::string_view tok : kGridOpTokens )
    {
        auto op = ParseGridOp( tok );
        ASSERT_TRUE( op.has_value() ) << tok;
        EXPECT_EQ( ToString( *op ), tok );
        EXPECT_TRUE( IsKnownGridOpToken( tok ) );
    }
}

TEST( GridOpCatalogTest, CaseInsensitiveParse )
{
    EXPECT_EQ( ParseGridOp( "calc_metrics" ), GridOp::CalcMetrics );
    EXPECT_EQ( ParseGridOp( "Fill_Wall_Struct" ), GridOp::FillWallStruct );
    EXPECT_EQ( ParseGridOp( "ALLOCATE_WALL_DIST" ), GridOp::AllocWallDist );
}

TEST( GridOpCatalogTest, UnknownTokenRejected )
{
    EXPECT_FALSE( ParseGridOp( "NOT_A_GRID_OP" ).has_value() );
    EXPECT_FALSE( IsKnownGridOpToken( "" ) );
    EXPECT_FALSE( IsKnownGridOpToken( "CALC_METRIC" ) ); // missing S
}

TEST( GridOpCatalogTest, TokenTableHasUniqueEntries )
{
    std::set< std::string > seen;
    for ( std::string_view tok : kGridOpTokens )
    {
        EXPECT_TRUE( seen.insert( std::string( tok ) ).second ) << tok;
    }
    EXPECT_EQ( seen.size(), kGridOpTokens.size() );
}

TEST( GridTypesTest, ParseGridFileTypeGridgen )
{
    EXPECT_EQ( ParseGridFileType( "gridgen" ), GridFileType::Gridgen );
    EXPECT_EQ( ParseGridFileType( "GridGen" ), GridFileType::Gridgen );
    EXPECT_EQ( ToString( GridFileType::Gridgen ), "gridgen" );
}

TEST( GridTypesTest, ParseGridTopology )
{
    EXPECT_EQ( ParseGridTopology( "u" ), GridTopology::Unstructured );
    EXPECT_EQ( ParseGridTopology( "UNSTRUCTURED" ), GridTopology::Unstructured );
    EXPECT_EQ( ParseGridTopology( "s" ), GridTopology::Structured );
    EXPECT_EQ( ParseGridTopology( "structured" ), GridTopology::Structured );
    EXPECT_EQ( ParseGridTopology( "other" ), GridTopology::Unknown );
}

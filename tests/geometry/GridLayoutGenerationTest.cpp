#include "DataBase.h"
#include "Dimension.h"
#include "GridMachine.h"
#include "Prj.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <system_error>
#include <vector>

namespace ONEFLOW
{

namespace
{

template < typename T >
T ReadBinary( std::ifstream & file )
{
    T value{};
    file.read( reinterpret_cast< char * >( & value ), sizeof( T ) );
    return value;
}

struct GridGenerationTestContext
{
    int dimension;
    std::string projectDir;
    std::filesystem::path outputDir;

    ~GridGenerationTestContext()
    {
        Dim::SetDimension( dimension );
        Prj::prjBaseDir = projectDir;
        std::error_code error;
        std::filesystem::remove_all( outputDir, error );
    }
};

TEST( GridLayoutGenerationTest, GeneratesMinimalRectangle )
{
    const std::filesystem::path layoutFile =
        std::filesystem::path( ONEFLOW_SOURCE_DIR ) /
        "tests/geometry/data/minimalRectangleLayout.txt";

    const std::filesystem::path outputDir =
        std::filesystem::temp_directory_path() / "oneflow-minimal-rectangle";
    std::filesystem::remove_all( outputDir );
    std::filesystem::create_directories( outputDir );

    const std::filesystem::path outputFile = outputDir / "minimalRectangle.ofl";
    const std::filesystem::path bcFile = outputDir / "minimalRectangle.inp";

    const int oldDimension = Dim::GetDimension();
    const std::string oldProjectDir = Prj::prjBaseDir;
    const GridGenerationTestContext restore{ oldDimension, oldProjectDir, outputDir };

    Dim::SetDimension( TWO_D );
    Prj::SetPrjBaseDir( outputDir.string() );

    SetDataString( "gridLayoutFileName", layoutFile.string() );
    SetDataString( "sourceGridFileName", outputFile.string() );
    SetDataString( "sourceGridBcName", bcFile.string() );
    SetDataString( "targetGridFileName", outputFile.string() );
    SetDataString( "part_uns_file", ( outputDir / "partition.uns" ).string() );

    GridMachine().Run( layoutFile.string() );

    ASSERT_TRUE( std::filesystem::is_regular_file( outputFile ) );

    std::ifstream file( outputFile, std::ios::binary );
    ASSERT_TRUE( file.is_open() );

    const int nZone = ReadBinary< int >( file );
    ASSERT_EQ( nZone, 1 );

    const int ni = ReadBinary< int >( file );
    const int nj = ReadBinary< int >( file );
    const int nk = ReadBinary< int >( file );

    EXPECT_EQ( ni, 3 );
    EXPECT_EQ( nj, 3 );
    EXPECT_EQ( nk, 1 );

    constexpr int nodeCount = 9;
    std::array< Real, nodeCount > x{};
    std::array< Real, nodeCount > y{};
    std::array< Real, nodeCount > z{};

    for ( Real & value : x ) value = ReadBinary< Real >( file );
    for ( Real & value : y ) value = ReadBinary< Real >( file );
    for ( Real & value : z ) value = ReadBinary< Real >( file );

    const std::array< Real, nodeCount > expectedX = {
        0.0, 0.5, 1.0,
        0.0, 0.5, 1.0,
        0.0, 0.5, 1.0
    };
    const std::array< Real, nodeCount > expectedY = {
        0.0, 0.0, 0.0,
        0.5, 0.5, 0.5,
        1.0, 1.0, 1.0
    };
    const std::array< Real, nodeCount > expectedZ = {};

    for ( int i = 0; i < nodeCount; ++ i )
    {
        EXPECT_DOUBLE_EQ( x[ i ], expectedX[ i ] );
        EXPECT_DOUBLE_EQ( y[ i ], expectedY[ i ] );
        EXPECT_DOUBLE_EQ( z[ i ], expectedZ[ i ] );
    }

    EXPECT_EQ( ( ni - 1 ) * ( nj - 1 ), 4 );

    file.close();
}


TEST( GridLayoutGenerationTest, GeneratesTwoConnectedThreeDimensionalCubes )
{
    const std::filesystem::path layoutFile =
        std::filesystem::path( ONEFLOW_SOURCE_DIR ) /
        "grid/two_connected_cubes/grid/twoConnectedCubesLayout.txt";

    const std::filesystem::path outputDir =
        std::filesystem::temp_directory_path() / "oneflow-two-connected-cubes-3d";
    std::filesystem::remove_all( outputDir );
    std::filesystem::create_directories( outputDir / "grid" );

    const std::filesystem::path outputFile = outputDir / "twoConnectedCubes.ofl";
    const std::filesystem::path bcFile = outputDir / "twoConnectedCubes.inp";

    const int oldDimension = Dim::GetDimension();
    const std::string oldProjectDir = Prj::prjBaseDir;
    const GridGenerationTestContext restore{ oldDimension, oldProjectDir, outputDir };

    Dim::SetDimension( THREE_D );
    Prj::SetPrjBaseDir( outputDir.string() );

    SetDataString( "gridLayoutFileName", layoutFile.string() );
    SetDataString( "sourceGridFileName", outputFile.string() );
    SetDataString( "sourceGridBcName", bcFile.string() );
    SetDataString( "targetGridFileName", outputFile.string() );
    SetDataString( "part_uns_file", ( outputDir / "partition.uns" ).string() );

    GridMachine().Run( layoutFile.string() );

    ASSERT_TRUE( std::filesystem::is_regular_file( outputFile ) );

    std::ifstream file( outputFile, std::ios::binary );
    ASSERT_TRUE( file.is_open() );

    ASSERT_EQ( ReadBinary< int >( file ), 2 );

    struct ZoneDimensions
    {
        int ni;
        int nj;
        int nk;
    };
    std::array< ZoneDimensions, 2 > dimensions{};
    for ( auto & zone : dimensions )
    {
        zone.ni = ReadBinary< int >( file );
        zone.nj = ReadBinary< int >( file );
        zone.nk = ReadBinary< int >( file );
        ASSERT_EQ( zone.ni, 2 );
        ASSERT_EQ( zone.nj, 2 );
        ASSERT_EQ( zone.nk, 2 );
        EXPECT_EQ( ( zone.ni - 1 ) * ( zone.nj - 1 ) * ( zone.nk - 1 ), 1 );
    }

    struct ZoneCoordinates
    {
        std::vector< Real > x;
        std::vector< Real > y;
        std::vector< Real > z;
    };
    std::array< ZoneCoordinates, 2 > zones;
    for ( std::size_t zoneIndex = 0; zoneIndex < zones.size(); ++ zoneIndex )
    {
        const auto & dims = dimensions[ zoneIndex ];
        const std::size_t nodeCount = static_cast< std::size_t >(
            dims.ni * dims.nj * dims.nk );
        zones[ zoneIndex ].x.resize( nodeCount );
        zones[ zoneIndex ].y.resize( nodeCount );
        zones[ zoneIndex ].z.resize( nodeCount );

        for ( Real & value : zones[ zoneIndex ].x ) value = ReadBinary< Real >( file );
        for ( Real & value : zones[ zoneIndex ].y ) value = ReadBinary< Real >( file );
        for ( Real & value : zones[ zoneIndex ].z ) value = ReadBinary< Real >( file );
    }

    for ( std::size_t zoneIndex = 0; zoneIndex < zones.size(); ++ zoneIndex )
    {
        const auto & zone = zones[ zoneIndex ];
        for ( const Real value : zone.x ) EXPECT_TRUE( std::isfinite( value ) );
        for ( const Real value : zone.y ) EXPECT_TRUE( std::isfinite( value ) );
        for ( const Real value : zone.z ) EXPECT_TRUE( std::isfinite( value ) );

        const Real expectedMinX = static_cast< Real >( zoneIndex );
        const Real expectedMaxX = expectedMinX + 1.0;

        EXPECT_DOUBLE_EQ(
            *std::min_element( zone.x.begin(), zone.x.end() ), expectedMinX );
        EXPECT_DOUBLE_EQ(
            *std::max_element( zone.x.begin(), zone.x.end() ), expectedMaxX );
        EXPECT_DOUBLE_EQ(
            *std::min_element( zone.y.begin(), zone.y.end() ), 0.0 );
        EXPECT_DOUBLE_EQ(
            *std::max_element( zone.y.begin(), zone.y.end() ), 1.0 );
        EXPECT_DOUBLE_EQ(
            *std::min_element( zone.z.begin(), zone.z.end() ), 0.0 );
        EXPECT_DOUBLE_EQ(
            *std::max_element( zone.z.begin(), zone.z.end() ), 1.0 );
    }

    int matchedInterfaceNodes = 0;
    for ( std::size_t i = 0; i < zones[ 0 ].x.size(); ++ i )
    {
        if ( zones[ 0 ].x[ i ] != 1.0 ) continue;

        bool matched = false;
        for ( std::size_t j = 0; j < zones[ 1 ].x.size(); ++ j )
        {
            if ( zones[ 1 ].x[ j ] == 1.0 &&
                 zones[ 1 ].y[ j ] == zones[ 0 ].y[ i ] &&
                 zones[ 1 ].z[ j ] == zones[ 0 ].z[ i ] )
            {
                matched = true;
                break;
            }
        }
        EXPECT_TRUE( matched );
        if ( matched ) ++ matchedInterfaceNodes;
    }
    EXPECT_EQ( matchedInterfaceNodes, 4 );

    file.close();
}

} // namespace

} // namespace ONEFLOW

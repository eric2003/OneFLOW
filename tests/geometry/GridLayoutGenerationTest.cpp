#include "DataBase.h"
#include "Dimension.h"
#include "GridMachine.h"
#include "Prj.h"

#include <gtest/gtest.h>

#include <array>
#include <filesystem>
#include <fstream>
#include <string>

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

    Dim::SetDimension( TWO_D );
    Prj::SetPrjBaseDir( outputDir.string() );

    SetDataString( "gridLayoutFileName", layoutFile.string() );
    SetDataString( "sourceGridFileName", outputFile.string() );
    SetDataString( "sourceGridBcName", bcFile.string() );
    SetDataString( "targetGridFileName", outputFile.string() );

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

    Dim::SetDimension( oldDimension );
    Prj::prjBaseDir = oldProjectDir;
    std::filesystem::remove_all( outputDir );
}

} // namespace

} // namespace ONEFLOW

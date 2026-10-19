#include "DataBase.h"
#include "Dimension.h"
#include "GridMachine.h"
#include "MDomain.h"
#include "MLine.h"
#include "Prj.h"
#include "SDomain.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
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

TEST( GridLayoutGenerationTest, GeneratesFlatPlateBoundaryLayerDistribution )
{
    const std::filesystem::path layoutFile =
        std::filesystem::path( ONEFLOW_SOURCE_DIR ) /
        "grid/laminar_flat_plate/grid/laminarFlatPlateLayout.txt";

    const std::filesystem::path outputDir =
        std::filesystem::temp_directory_path() / "oneflow-flat-plate-boundary-layer";
    std::filesystem::remove_all( outputDir );
    std::filesystem::create_directories( outputDir );

    const std::filesystem::path outputFile =
        outputDir / "laminarFlatPlate.generated.grd";
    const std::filesystem::path bcFile =
        outputDir / "laminarFlatPlate.inp";

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

    ASSERT_EQ( ReadBinary< int >( file ), 2 );

    struct ZoneDimensions
    {
        int ni;
        int nj;
        int nk;
    };
    const std::array< ZoneDimensions, 2 > expectedDimensions = {{
        { 33, 85, 1 },
        { 41, 85, 1 }
    }};
    std::array< ZoneDimensions, 2 > dimensions{};
    for ( std::size_t zoneIndex = 0; zoneIndex < dimensions.size(); ++ zoneIndex )
    {
        auto & zone = dimensions[ zoneIndex ];
        zone.ni = ReadBinary< int >( file );
        zone.nj = ReadBinary< int >( file );
        zone.nk = ReadBinary< int >( file );

        EXPECT_EQ( zone.ni, expectedDimensions[ zoneIndex ].ni );
        EXPECT_EQ( zone.nj, expectedDimensions[ zoneIndex ].nj );
        EXPECT_EQ( zone.nk, expectedDimensions[ zoneIndex ].nk );
    }

    constexpr int nearWallPointCount = 41;
    constexpr Real nearWallDistance = 0.02;
    constexpr Real domainHeight = 4.0;

    struct ZoneCoordinates
    {
        std::vector< Real > x;
        std::vector< Real > y;
        std::vector< Real > z;
    };
    std::array< ZoneCoordinates, 2 > zones;
    for ( std::size_t zoneIndex = 0; zoneIndex < zones.size(); ++ zoneIndex )
    {
        const std::size_t nodeCount = static_cast< std::size_t >(
            dimensions[ zoneIndex ].ni * dimensions[ zoneIndex ].nj );
        zones[ zoneIndex ].x.resize( nodeCount );
        zones[ zoneIndex ].y.resize( nodeCount );
        zones[ zoneIndex ].z.resize( nodeCount );

        for ( Real & value : zones[ zoneIndex ].x ) value = ReadBinary< Real >( file );
        for ( Real & value : zones[ zoneIndex ].y ) value = ReadBinary< Real >( file );
        for ( Real & value : zones[ zoneIndex ].z ) value = ReadBinary< Real >( file );
    }

    for ( const auto & zone : zones )
    {
        for ( const Real value : zone.x ) EXPECT_TRUE( std::isfinite( value ) );
        for ( const Real value : zone.y ) EXPECT_TRUE( std::isfinite( value ) );
        for ( const Real value : zone.z ) EXPECT_TRUE( std::isfinite( value ) );
    }

    const auto nodeIndex = []( int i, int j, int ni )
    {
        return static_cast< std::size_t >( j * ni + i );
    };
    for ( int j = 0; j < dimensions[ 0 ].nj; ++ j )
    {
        const std::size_t leftInterfaceNode =
            nodeIndex( dimensions[ 0 ].ni - 1, j, dimensions[ 0 ].ni );
        const std::size_t rightInterfaceNode =
            nodeIndex( 0, j, dimensions[ 1 ].ni );
        EXPECT_NEAR( zones[ 0 ].x[ leftInterfaceNode ], 0.0, 1.0e-12 );
        EXPECT_NEAR( zones[ 1 ].x[ rightInterfaceNode ], 0.0, 1.0e-12 );
        EXPECT_NEAR(
            zones[ 0 ].y[ leftInterfaceNode ],
            zones[ 1 ].y[ rightInterfaceNode ],
            1.0e-12 );
        EXPECT_NEAR(
            zones[ 0 ].z[ leftInterfaceNode ],
            zones[ 1 ].z[ rightInterfaceNode ],
            1.0e-12 );
    }

    for ( std::size_t zoneIndex = 0; zoneIndex < zones.size(); ++ zoneIndex )
    {
        const auto & zone = zones[ zoneIndex ];
        const int ni = dimensions[ zoneIndex ].ni;
        const int nj = dimensions[ zoneIndex ].nj;
        Real totalCellArea = 0.0;
        for ( int j = 0; j < nj - 1; ++ j )
        {
            for ( int i = 0; i < ni - 1; ++ i )
            {
                const std::size_t lowerLeft = nodeIndex( i, j, ni );
                const std::size_t lowerRight = nodeIndex( i + 1, j, ni );
                const std::size_t upperRight = nodeIndex( i + 1, j + 1, ni );
                const std::size_t upperLeft = nodeIndex( i, j + 1, ni );
                const Real cellArea = 0.5 * (
                    zone.x[ lowerLeft ] * zone.y[ lowerRight ] -
                    zone.y[ lowerLeft ] * zone.x[ lowerRight ] +
                    zone.x[ lowerRight ] * zone.y[ upperRight ] -
                    zone.y[ lowerRight ] * zone.x[ upperRight ] +
                    zone.x[ upperRight ] * zone.y[ upperLeft ] -
                    zone.y[ upperRight ] * zone.x[ upperLeft ] +
                    zone.x[ upperLeft ] * zone.y[ lowerLeft ] -
                    zone.y[ upperLeft ] * zone.x[ lowerLeft ] );
                EXPECT_GT( cellArea, 0.0 );
                totalCellArea += cellArea;
            }
        }
        EXPECT_NEAR( totalCellArea, 8.0, 1.0e-10 );
    }

    const Real uniformStep =
        nearWallDistance / static_cast< Real >( nearWallPointCount - 1 );
    for ( std::size_t zoneIndex = 0; zoneIndex < zones.size(); ++ zoneIndex )
    {
        const auto & zone = zones[ zoneIndex ];
        const int ni = dimensions[ zoneIndex ].ni;
        const int nj = dimensions[ zoneIndex ].nj;

        for ( int i = 0; i < ni; ++ i )
        {
            const auto yAt = [ &zone, ni, i ]( int j ) -> Real
            {
                return zone.y[ static_cast< std::size_t >( j * ni + i ) ];
            };

            const int wallIndex =
                std::abs( yAt( 0 ) ) <= std::abs( yAt( nj - 1 ) ) ? 0 : nj - 1;
            const int outwardDirection = wallIndex == 0 ? 1 : -1;
            const auto yFromWall = [ &yAt, wallIndex, outwardDirection ]( int layer ) -> Real
            {
                return yAt( wallIndex + outwardDirection * layer );
            };

            for ( int layer = 0; layer < nearWallPointCount; ++ layer )
            {
                const Real expectedY = uniformStep * static_cast< Real >( layer );
                EXPECT_NEAR( yFromWall( layer ), expectedY, 1.0e-12 );
            }

            Real previousStep = yFromWall( 1 ) - yFromWall( 0 );
            ASSERT_GT( previousStep, 0.0 );
            EXPECT_NEAR( previousStep, uniformStep, 1.0e-12 );

            const Real firstOuterStep =
                yFromWall( nearWallPointCount ) - yFromWall( nearWallPointCount - 1 );
            EXPECT_NEAR( firstOuterStep, previousStep, 1.0e-12 );

            const Real secondOuterStep =
                yFromWall( nearWallPointCount + 1 ) - yFromWall( nearWallPointCount );
            EXPECT_GT( secondOuterStep, firstOuterStep );
            const Real outerRatio = secondOuterStep / firstOuterStep;
            EXPECT_GT( outerRatio, 1.0 );
            previousStep = secondOuterStep;

            for ( int layer = nearWallPointCount + 2; layer < nj; ++ layer )
            {
                const Real step = yFromWall( layer ) - yFromWall( layer - 1 );
                EXPECT_GT( step, previousStep );
                EXPECT_NEAR( step / previousStep, outerRatio, 1.0e-10 );
                previousStep = step;
            }

            EXPECT_NEAR( yFromWall( nj - 1 ), domainHeight, 1.0e-12 );
            EXPECT_NEAR( yFromWall( 0 ), 0.0, 1.0e-12 );
        }
    }

    file.close();
}

TEST( GridLayoutGenerationTest, AssemblesRectangleFromTwoLocallyGeneratedBlocks )
{
    const std::filesystem::path layoutFile =
        std::filesystem::path( ONEFLOW_SOURCE_DIR ) /
        "grid/composite_boundary/grid.txt";
    const std::filesystem::path outputDir =
        std::filesystem::temp_directory_path() / "oneflow-composite-boundary";
    std::filesystem::remove_all( outputDir );
    std::filesystem::create_directories( outputDir );

    const std::filesystem::path outputFile = outputDir / "compositeBoundary.grd";
    const std::filesystem::path bcFile = outputDir / "compositeBoundary.inp";
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
    ASSERT_EQ( ReadBinary< int >( file ), 2 );

    struct Zone
    {
        int ni = 0;
        int nj = 0;
        int nk = 0;
        std::vector< Real > x;
        std::vector< Real > y;
        std::vector< Real > z;
    };
    std::array< Zone, 2 > zones{};
    for ( Zone & zone : zones )
    {
        zone.ni = ReadBinary< int >( file );
        zone.nj = ReadBinary< int >( file );
        zone.nk = ReadBinary< int >( file );
        ASSERT_EQ( zone.ni * zone.nj, 4 );
        EXPECT_EQ( zone.nk, 1 );
    }

    for ( Zone & zone : zones )
    {
        const std::size_t nodeCount =
            static_cast< std::size_t >( zone.ni * zone.nj );
        zone.x.resize( nodeCount );
        zone.y.resize( nodeCount );
        zone.z.resize( nodeCount );
        for ( Real & value : zone.x ) value = ReadBinary< Real >( file );
        for ( Real & value : zone.y ) value = ReadBinary< Real >( file );
        for ( Real & value : zone.z ) value = ReadBinary< Real >( file );
    }

    struct Point2D
    {
        Real x;
        Real y;
    };
    const auto containsPoint = []( const Zone & zone, Real x, Real y )
    {
        for ( std::size_t node = 0; node < zone.x.size(); ++ node )
        {
            if ( std::abs( zone.x[ node ] - x ) <= 1.0e-12 &&
                 std::abs( zone.y[ node ] - y ) <= 1.0e-12 )
            {
                return true;
            }
        }
        return false;
    };

    const std::array< std::array< Point2D, 4 >, 2 > expectedCorners = {{
        {{{ 0.0, 0.0 }, { 4.0, 0.0 }, { 4.0, 1.0 }, { 0.0, 1.0 }}},
        {{{ 0.0, 1.0 }, { 4.0, 1.0 }, { 4.0, 2.0 }, { 0.0, 2.0 }}}
    }};
    Real totalArea = 0.0;
    for ( std::size_t zoneIndex = 0; zoneIndex < zones.size(); ++ zoneIndex )
    {
        const Zone & zone = zones[ zoneIndex ];
        for ( std::size_t node = 0; node < zone.x.size(); ++ node )
        {
            EXPECT_TRUE( std::isfinite( zone.x[ node ] ) );
            EXPECT_TRUE( std::isfinite( zone.y[ node ] ) );
            EXPECT_DOUBLE_EQ( zone.z[ node ], 0.0 );
            EXPECT_GE( zone.x[ node ], 0.0 );
            EXPECT_LE( zone.x[ node ], 4.0 );
            EXPECT_GE( zone.y[ node ], 0.0 );
            EXPECT_LE( zone.y[ node ], 2.0 );
        }

        for ( const Point2D & corner : expectedCorners[ zoneIndex ] )
        {
            EXPECT_TRUE( containsPoint( zone, corner.x, corner.y ) );
        }

        Real twiceSignedArea = 0.0;
        for ( std::size_t node = 0; node < zone.x.size(); ++ node )
        {
            const std::size_t next = ( node + 1 ) % zone.x.size();
            twiceSignedArea +=
                zone.x[ node ] * zone.y[ next ] -
                zone.y[ node ] * zone.x[ next ];
        }
        const Real area = 0.5 * std::abs( twiceSignedArea );
        EXPECT_NEAR( area, 4.0, 1.0e-12 );
        totalArea += area;
    }
    EXPECT_NEAR( totalArea, 8.0, 1.0e-12 );

    int sharedInterfaceNodeCount = 0;
    for ( std::size_t firstNode = 0; firstNode < zones[ 0 ].x.size(); ++ firstNode )
    {
        for ( std::size_t secondNode = 0; secondNode < zones[ 1 ].x.size(); ++ secondNode )
        {
            if ( std::abs( zones[ 0 ].x[ firstNode ] - zones[ 1 ].x[ secondNode ] ) <= 1.0e-12 &&
                 std::abs( zones[ 0 ].y[ firstNode ] - 1.0 ) <= 1.0e-12 &&
                 std::abs( zones[ 1 ].y[ secondNode ] - 1.0 ) <= 1.0e-12 )
            {
                EXPECT_NEAR( zones[ 0 ].x[ firstNode ], zones[ 1 ].x[ secondNode ], 1.0e-12 );
                EXPECT_NEAR( zones[ 0 ].z[ firstNode ], zones[ 1 ].z[ secondNode ], 1.0e-12 );
                ++ sharedInterfaceNodeCount;
            }
        }
    }
    EXPECT_EQ( sharedInterfaceNodeCount, 2 );

    file.close();
}

TEST( DomDataTest, RebuildsControlPointCandidatesFromCurrentDomainIncidence )
{
    DomData domainData;
    domainData.pointToDomainMap[ 101 ].insert( 1 );
    domainData.pointToDomainMap[ 205 ].insert( 1 );
    domainData.pointToDomainMap[ 205 ].insert( 2 );

    domainData.ConstructCtrlPoint();

    EXPECT_EQ( domainData.candidate_ctrlpoints, ( IntField{ 101 } ) );

    domainData.pointToDomainMap.clear();
    domainData.pointToDomainMap[ 309 ].insert( 3 );
    domainData.ConstructCtrlPoint();

    EXPECT_EQ( domainData.candidate_ctrlpoints, ( IntField{ 309 } ) );
}

TEST( DomDataTest, RebuildsBoundaryPointCandidatesFromCurrentBoundaryLines )
{
    DomData domainData;
    domainData.candidate_bcpoints = { 101, 205 };
    domainData.lineToDomainMap.clear();

    domainData.ConstructBcPoint();

    EXPECT_TRUE( domainData.candidate_bcpoints.empty() );
}

TEST( MLineTest, RebuildsEndpointCandidatesAfterTopologyChanges )
{
    CoorMap coordinates;
    MLine line( &coordinates );
    line.lineList = { 10, 20 };

    LinkField initialPointIds;
    initialPointIds.push_back( IntField{ 101, 205 } );
    initialPointIds.push_back( IntField{ 205, 309 } );
    line.ConstructPointToLineMap( initialPointIds );
    line.ConstructEndpointCandidates();

    EXPECT_EQ( line.endpointCandidates, ( IntField{ 101, 309 } ) );

    LinkField updatedPointIds;
    updatedPointIds.push_back( IntField{ 101, 205 } );
    updatedPointIds.push_back( IntField{ 205, 407 } );
    line.ConstructPointToLineMap( updatedPointIds );
    line.ConstructEndpointCandidates();

    EXPECT_EQ( line.endpointCandidates, ( IntField{ 101, 407 } ) );
    EXPECT_EQ( line.pointToLine.count( 309 ), 0u );
    EXPECT_EQ( line.pointToLine.at( 407 ).size(), 1u );
}

TEST( SDomainTest, RebuildsOwnedIncidenceMapsWithoutStaleEntries )
{
    CoorMap coordinates;
    SDomain domain( &coordinates );
    domain.pointToPointMap[ 101 ].insert( 205 );
    domain.pointToDomainMap[ 101 ].insert( 7 );
    domain.lineToDomainMap[ 10 ].insert( 7 );

    domain.ConstructPointToPointMap();
    domain.ConstructPointToDomainMap();
    domain.ConstructLineToDomainMap();

    EXPECT_TRUE( domain.pointToPointMap.empty() );
    EXPECT_TRUE( domain.pointToDomainMap.empty() );
    EXPECT_TRUE( domain.lineToDomainMap.empty() );
}

TEST( MLineTest, RebuildsOwnedIncidenceMapsWithoutStaleEntries )
{
    CoorMap coordinates;
    MLine line( &coordinates );
    line.pointToPointMap[ 101 ].insert( 205 );
    line.lineToDomainMap[ 10 ].insert( 7 );

    line.ConstructPointToPointMap();
    line.ConstructLineToDomainMap();

    EXPECT_TRUE( line.pointToPointMap.empty() );
    EXPECT_TRUE( line.lineToDomainMap.empty() );
}

TEST( SDomainTest, FindsUniqueSharedEndpointBetweenAdjacentEdges )
{
    CoorMap coordinates;
    SDomain domain( &coordinates );
    MLine firstEdge( &coordinates );
    MLine secondEdge( &coordinates );
    firstEdge.endpointCandidates = { 101, 205 };
    secondEdge.endpointCandidates = { 205, 309 };

    int commonPoint = -1;
    domain.GetCommonPoint( &firstEdge, &secondEdge, commonPoint );

    EXPECT_EQ( commonPoint, 205 );
}

TEST( SDomainTest, RejectsAdjacentEdgesWithMultipleSharedEndpointCandidates )
{
    CoorMap coordinates;
    SDomain domain( &coordinates );
    MLine firstEdge( &coordinates );
    MLine secondEdge( &coordinates );
    firstEdge.endpointCandidates = { 101, 205 };
    secondEdge.endpointCandidates = { 101, 205 };

    int commonPoint = -1;
    EXPECT_THROW(
        domain.GetCommonPoint( &firstEdge, &secondEdge, commonPoint ),
        std::runtime_error );
}

TEST( SDomainTest, RejectsAdjacentEdgesWithoutSharedEndpointCandidates )
{
    CoorMap coordinates;
    SDomain domain( &coordinates );
    MLine firstEdge( &coordinates );
    MLine secondEdge( &coordinates );
    firstEdge.endpointCandidates = { 101, 205 };
    secondEdge.endpointCandidates = { 309, 407 };

    int commonPoint = -1;
    EXPECT_THROW(
        domain.GetCommonPoint( &firstEdge, &secondEdge, commonPoint ),
        std::runtime_error );
}

TEST( MDomainTest, RebuildsAggregateIncidenceMapsWithoutStaleEntries )
{
    MDomain domain;
    domain.lineToDomainMap[ 10 ].insert( 7 );
    domain.pointToDomainMap[ 101 ].insert( 7 );
    domain.pointToPointMap[ 101 ].insert( 205 );
    domain.pointToLineMap[ 101 ].insert( 10 );

    domain.ConstructMultiLineToDomainMap();
    domain.ConstructMultiPointToDomainMap();
    domain.ConstructMultiPointToPointMap();
    domain.ConstructPointToLineMap();

    EXPECT_TRUE( domain.lineToDomainMap.empty() );
    EXPECT_TRUE( domain.pointToDomainMap.empty() );
    EXPECT_TRUE( domain.pointToPointMap.empty() );
    EXPECT_TRUE( domain.pointToLineMap.empty() );
}

TEST( MDomainTest, ReportsStalledControlPointResolution )
{
    CoorMap coordinates;

    CalcCoor lowerLeft;
    lowerLeft.SetCoor( 10, 20, 0 );
    coordinates.emplace( 101, lowerLeft );

    CalcCoor lowerRight;
    lowerRight.SetCoor( 20, 20, 0 );
    coordinates.emplace( 205, lowerRight );

    MDomain domain;
    domain.coorMap = &coordinates;

    auto subdomain = std::make_unique< SDomain >( &coordinates );
    subdomain->ctrlpoints.push_back( 101 );
    subdomain->ctrlpoints.push_back( 205 );
    subdomain->ctrlpoints.push_back( 309 );
    subdomain->ctrlpoints.push_back( 407 );
    domain.sDomainList.push_back( std::move( subdomain ) );

    EXPECT_THROW(
        domain.CalcSubDomainCtrlCoor(),
        std::runtime_error );
}

TEST( MLineTest, BuildsPointToLineIncidenceForCompositeLine )
{
    CoorMap coordinates;
    MLine line( &coordinates );
    line.lineList = { 10, 20 };

    LinkField pointIdLink;
    pointIdLink.push_back( IntField{ 101, 205 } );
    pointIdLink.push_back( IntField{ 205, 309 } );

    line.ConstructPointToLineMap( pointIdLink );
    line.ConstructEndpointCandidates();

    EXPECT_EQ( line.endpointCandidates, ( IntField{ 101, 309 } ) );
    ASSERT_EQ( line.pointToLine.size(), 3u );
    ASSERT_EQ( line.pointToLine.at( 101 ).size(), 1u );
    EXPECT_TRUE( line.pointToLine.at( 101 ).count( 10 ) );
    ASSERT_EQ( line.pointToLine.at( 205 ).size(), 2u );
    EXPECT_TRUE( line.pointToLine.at( 205 ).count( 10 ) );
    EXPECT_TRUE( line.pointToLine.at( 205 ).count( 20 ) );
    ASSERT_EQ( line.pointToLine.at( 309 ).size(), 1u );
    EXPECT_TRUE( line.pointToLine.at( 309 ).count( 20 ) );
}

TEST( SDomainTest, ReportsBlockedWhenFewerThanThreeControlPointsAreKnown )
{
    CoorMap coordinates;

    CalcCoor lowerLeft;
    lowerLeft.SetCoor( 10, 20, 0 );
    coordinates.emplace( 101, lowerLeft );

    CalcCoor lowerRight;
    lowerRight.SetCoor( 20, 20, 0 );
    coordinates.emplace( 205, lowerRight );

    SDomain domain( &coordinates );
    domain.ctrlpoints.push_back( 101 );
    domain.ctrlpoints.push_back( 205 );
    domain.ctrlpoints.push_back( 309 );
    domain.ctrlpoints.push_back( 407 );

    EXPECT_EQ(
        domain.CalcSingleDomainCoor(),
        ControlPointResolution::Blocked );
}

TEST( SDomainTest, ReportsProgressWhenInferringMissingControlPoint )
{
    CoorMap coordinates;

    CalcCoor lowerLeft;
    lowerLeft.SetCoor( 10, 20, 0 );
    coordinates.emplace( 101, lowerLeft );

    CalcCoor lowerRight;
    lowerRight.SetCoor( 20, 20, 0 );
    coordinates.emplace( 205, lowerRight );

    CalcCoor upperLeft;
    upperLeft.SetCoor( 10, 30, 0 );
    coordinates.emplace( 407, upperLeft );

    SDomain domain( &coordinates );
    domain.ctrlpoints.push_back( 101 );
    domain.ctrlpoints.push_back( 205 );
    domain.ctrlpoints.push_back( 309 );
    domain.ctrlpoints.push_back( 407 );

    EXPECT_EQ(
        domain.CalcSingleDomainCoor(),
        ControlPointResolution::Progress );

    const auto inferredCorner = coordinates.find( 309 );
    ASSERT_NE( inferredCorner, coordinates.end() );
    EXPECT_EQ( inferredCorner->second.i, 20 );
    EXPECT_EQ( inferredCorner->second.j, 30 );
    EXPECT_EQ( inferredCorner->second.k, 0 );
}

TEST( SDomainTest, CompletesMissingControlPointUsingGlobalPointIds )
{
    CoorMap coordinates;

    CalcCoor lowerLeft;
    lowerLeft.SetCoor( 10, 20, 0 );
    coordinates.emplace( 101, lowerLeft );

    CalcCoor lowerRight;
    lowerRight.SetCoor( 20, 20, 0 );
    coordinates.emplace( 205, lowerRight );

    CalcCoor upperLeft;
    upperLeft.SetCoor( 10, 30, 0 );
    coordinates.emplace( 407, upperLeft );

    SDomain domain( &coordinates );
    domain.ctrlpoints.push_back( 101 );
    domain.ctrlpoints.push_back( 205 );
    domain.ctrlpoints.push_back( 309 );
    domain.ctrlpoints.push_back( 407 );

    IntField knownCorners;
    knownCorners.push_back( 1 );
    knownCorners.push_back( 1 );
    knownCorners.push_back( 0 );
    knownCorners.push_back( 1 );

    domain.SetRemainingCtrlPoint( knownCorners );

    const auto missingCorner = coordinates.find( 309 );
    ASSERT_NE( missingCorner, coordinates.end() );
    EXPECT_EQ( missingCorner->second.i, 20 );
    EXPECT_EQ( missingCorner->second.j, 30 );
    EXPECT_EQ( missingCorner->second.k, 0 );
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



TEST( GridLayoutGenerationTest, GeneratesRectangleWithCompositeBoundaryEdge )
{
    const std::filesystem::path outputDir =
        std::filesystem::temp_directory_path() / "oneflow-composite-boundary-edge";
    std::filesystem::remove_all( outputDir );
    std::filesystem::create_directories( outputDir );

    const std::filesystem::path layoutFile = outputDir / "compositeRectangleLayout.txt";
    const std::filesystem::path outputFile = outputDir / "compositeRectangle.ofl";
    const std::filesystem::path bcFile = outputDir / "compositeRectangle.inp";

    std::ofstream layout( layoutFile );
    ASSERT_TRUE( layout.is_open() );
    layout <<
        "Point(1)={0,0,0}\n"
        "Point(2)={1,0,0}\n"
        "Point(3)={1,1,0}\n"
        "Point(4)={0.5,1,0}\n"
        "Point(5)={0,1,0}\n"
        "Line(1)={1,2}\n"
        "Line(2)={2,3}\n"
        "Line(3)={3,4}\n"
        "Line(4)={4,5}\n"
        "Line(5)={5,1}\n"
        "Dim(1)={3}\n"
        "Dim(2)={3}\n"
        "Dim(3)={2}\n"
        "Dim(4)={2}\n"
        "Dim(5)={3}\n"
        "Ds(1)={d,-1,-1}\n"
        "Ds(2)={d,-1,-1}\n"
        "Ds(3)={d,-1,-1}\n"
        "Ds(4)={d,-1,-1}\n"
        "Ds(5)={d,-1,-1}\n"
        "Boundary(1)={1}\n"
        "Boundary(2)={2}\n"
        "Boundary(3)={3}\n"
        "Boundary(4)={3}\n"
        "Boundary(5)={4}\n"
        "Add L2F 1 1 1\n"
        "Add L2F 1 2 2\n"
        "Add L2F 1 3 3\n"
        "Add L2F 1 3 4\n"
        "Add L2F 1 4 5\n"
        "Add F2B 1 1 1\n";
    layout.close();

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
    ASSERT_EQ( ReadBinary< int >( file ), 1 );
    EXPECT_EQ( ReadBinary< int >( file ), 3 );
    EXPECT_EQ( ReadBinary< int >( file ), 3 );
    EXPECT_EQ( ReadBinary< int >( file ), 1 );

    constexpr int nodeCount = 9;
    std::array< Real, nodeCount > x{};
    std::array< Real, nodeCount > y{};
    std::array< Real, nodeCount > z{};
    for ( Real & value : x ) value = ReadBinary< Real >( file );
    for ( Real & value : y ) value = ReadBinary< Real >( file );
    for ( Real & value : z ) value = ReadBinary< Real >( file );

    for ( int i = 0; i < nodeCount; ++ i )
    {
        EXPECT_TRUE( std::isfinite( x[ i ] ) );
        EXPECT_TRUE( std::isfinite( y[ i ] ) );
        EXPECT_TRUE( std::isfinite( z[ i ] ) );
    }

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

    bool foundCompositeVertex = false;
    for ( int i = 0; i < nodeCount; ++ i )
    {
        if ( std::abs( x[ i ] - 0.5 ) < 1.0e-12 &&
             std::abs( y[ i ] - 1.0 ) < 1.0e-12 )
        {
            foundCompositeVertex = true;
            break;
        }
    }
    EXPECT_TRUE( foundCompositeVertex );

    Real totalCellArea = 0.0;
    for ( int j = 0; j < 2; ++ j )
    {
        for ( int i = 0; i < 2; ++ i )
        {
            const int lowerLeft = j * 3 + i;
            const int lowerRight = lowerLeft + 1;
            const int upperLeft = lowerLeft + 3;
            const int upperRight = upperLeft + 1;

            const Real cellArea =
                0.5 * (
                    x[ lowerLeft ] * y[ lowerRight ] -
                    y[ lowerLeft ] * x[ lowerRight ] +
                    x[ lowerRight ] * y[ upperRight ] -
                    y[ lowerRight ] * x[ upperRight ] +
                    x[ upperRight ] * y[ upperLeft ] -
                    y[ upperRight ] * x[ upperLeft ] +
                    x[ upperLeft ] * y[ lowerLeft ] -
                    y[ upperLeft ] * x[ lowerLeft ] );

            EXPECT_GT( cellArea, 0.0 );
            EXPECT_DOUBLE_EQ( cellArea, 0.25 );
            totalCellArea += cellArea;
        }
    }

    EXPECT_DOUBLE_EQ( totalCellArea, 1.0 );
    file.close();
}

TEST( GridLayoutGenerationTest, GeneratesRectangleWithTwoCompositeBoundaryEdges )
{
    const std::filesystem::path outputDir =
        std::filesystem::temp_directory_path() / "oneflow-two-composite-boundary-edges";
    std::filesystem::remove_all( outputDir );
    std::filesystem::create_directories( outputDir );

    const std::filesystem::path layoutFile = outputDir / "twoCompositeRectangleLayout.txt";
    const std::filesystem::path outputFile = outputDir / "twoCompositeRectangle.ofl";
    const std::filesystem::path bcFile = outputDir / "twoCompositeRectangle.inp";

    std::ofstream layout( layoutFile );
    ASSERT_TRUE( layout.is_open() );
    layout <<
        "Point(1)={0,0,0}\n"
        "Point(2)={1,0,0}\n"
        "Point(3)={1,1,0}\n"
        "Point(4)={0.5,1,0}\n"
        "Point(5)={0,1,0}\n"
        "Point(6)={0.5,0,0}\n"
        "Line(1)={6,1}\n"
        "Line(2)={3,2}\n"
        "Line(3)={3,4}\n"
        "Line(4)={5,4}\n"
        "Line(5)={1,5}\n"
        "Line(6)={6,2}\n"
        "Dim(1)={2}\n"
        "Dim(2)={3}\n"
        "Dim(3)={2}\n"
        "Dim(4)={2}\n"
        "Dim(5)={3}\n"
        "Dim(6)={2}\n"
        "Ds(1)={d,-1,-1}\n"
        "Ds(2)={d,-1,-1}\n"
        "Ds(3)={d,-1,-1}\n"
        "Ds(4)={d,-1,-1}\n"
        "Ds(5)={d,-1,-1}\n"
        "Ds(6)={d,-1,-1}\n"
        "Boundary(1)={1}\n"
        "Boundary(2)={2}\n"
        "Boundary(3)={3}\n"
        "Boundary(4)={3}\n"
        "Boundary(5)={4}\n"
        "Boundary(6)={1}\n"
        "Add L2F 1 1 1\n"
        "Add L2F 1 2 2\n"
        "Add L2F 1 3 3\n"
        "Add L2F 1 3 4\n"
        "Add L2F 1 4 5\n"
        "Add L2F 1 1 6\n"
        "Add F2B 1 1 1\n";
    layout.close();

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
    ASSERT_EQ( ReadBinary< int >( file ), 1 );
    EXPECT_EQ( ReadBinary< int >( file ), 3 );
    EXPECT_EQ( ReadBinary< int >( file ), 3 );
    EXPECT_EQ( ReadBinary< int >( file ), 1 );

    constexpr int nodeCount = 9;
    std::array< Real, nodeCount > x{};
    std::array< Real, nodeCount > y{};
    std::array< Real, nodeCount > z{};
    for ( Real & value : x ) value = ReadBinary< Real >( file );
    for ( Real & value : y ) value = ReadBinary< Real >( file );
    for ( Real & value : z ) value = ReadBinary< Real >( file );

    for ( int i = 0; i < nodeCount; ++ i )
    {
        EXPECT_TRUE( std::isfinite( x[ i ] ) );
        EXPECT_TRUE( std::isfinite( y[ i ] ) );
        EXPECT_TRUE( std::isfinite( z[ i ] ) );
    }

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

    bool foundCompositeVertex = false;
    for ( int i = 0; i < nodeCount; ++ i )
    {
        if ( std::abs( x[ i ] - 0.5 ) < 1.0e-12 &&
             std::abs( y[ i ] - 1.0 ) < 1.0e-12 )
        {
            foundCompositeVertex = true;
            break;
        }
    }
    EXPECT_TRUE( foundCompositeVertex );

    Real totalCellArea = 0.0;
    for ( int j = 0; j < 2; ++ j )
    {
        for ( int i = 0; i < 2; ++ i )
        {
            const int lowerLeft = j * 3 + i;
            const int lowerRight = lowerLeft + 1;
            const int upperLeft = lowerLeft + 3;
            const int upperRight = upperLeft + 1;

            const Real cellArea =
                0.5 * (
                    x[ lowerLeft ] * y[ lowerRight ] -
                    y[ lowerLeft ] * x[ lowerRight ] +
                    x[ lowerRight ] * y[ upperRight ] -
                    y[ lowerRight ] * x[ upperRight ] +
                    x[ upperRight ] * y[ upperLeft ] -
                    y[ upperRight ] * x[ upperLeft ] +
                    x[ upperLeft ] * y[ lowerLeft ] -
                    y[ upperLeft ] * x[ lowerLeft ] );

            EXPECT_GT( cellArea, 0.0 );
            EXPECT_DOUBLE_EQ( cellArea, 0.25 );
            totalCellArea += cellArea;
        }
    }

    EXPECT_DOUBLE_EQ( totalCellArea, 1.0 );
    file.close();
}

TEST( GridLayoutGenerationTest, ReordersThreeCompositeBoundarySegmentsRegardlessOfDirection )
{
    const std::filesystem::path outputDir =
        std::filesystem::temp_directory_path() / "oneflow-three-composite-boundary-segments";
    std::filesystem::remove_all( outputDir );
    std::filesystem::create_directories( outputDir );

    const std::filesystem::path layoutFile = outputDir / "threeSegmentRectangleLayout.txt";
    const std::filesystem::path outputFile = outputDir / "threeSegmentRectangle.ofl";
    const std::filesystem::path bcFile = outputDir / "threeSegmentRectangle.inp";

    std::ofstream layout( layoutFile );
    ASSERT_TRUE( layout.is_open() );
    layout <<
        "Point(1)={0,0,0}\n"
        "Point(2)={1,0,0}\n"
        "Point(3)={1,1,0}\n"
        "Point(4)={0,1,0}\n"
        "Point(5)={0.3333333333333333,0,0}\n"
        "Point(6)={0.6666666666666666,0,0}\n"
        "Line(1)={1,5}\n"
        "Line(2)={2,3}\n"
        "Line(3)={3,4}\n"
        "Line(4)={6,2}\n"
        "Line(5)={4,1}\n"
        "Line(6)={6,5}\n"
        "Dim(1)={2}\n"
        "Dim(2)={4}\n"
        "Dim(3)={4}\n"
        "Dim(4)={2}\n"
        "Dim(5)={4}\n"
        "Dim(6)={2}\n"
        "Ds(1)={d,-1,-1}\n"
        "Ds(2)={d,-1,-1}\n"
        "Ds(3)={d,-1,-1}\n"
        "Ds(4)={d,-1,-1}\n"
        "Ds(5)={d,-1,-1}\n"
        "Ds(6)={d,-1,-1}\n"
        "Boundary(1)={1}\n"
        "Boundary(2)={2}\n"
        "Boundary(3)={3}\n"
        "Boundary(4)={1}\n"
        "Boundary(5)={4}\n"
        "Boundary(6)={1}\n"
        "Add L2F 1 1 1\n"
        "Add L2F 1 2 2\n"
        "Add L2F 1 3 3\n"
        "Add L2F 1 1 4\n"
        "Add L2F 1 4 5\n"
        "Add L2F 1 1 6\n"
        "Add F2B 1 1 1\n";
    layout.close();

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
    ASSERT_EQ( ReadBinary< int >( file ), 1 );
    EXPECT_EQ( ReadBinary< int >( file ), 4 );
    EXPECT_EQ( ReadBinary< int >( file ), 4 );
    EXPECT_EQ( ReadBinary< int >( file ), 1 );

    constexpr int nodeCount = 16;
    std::array< Real, nodeCount > x{};
    std::array< Real, nodeCount > y{};
    std::array< Real, nodeCount > z{};
    for ( Real & value : x ) value = ReadBinary< Real >( file );
    for ( Real & value : y ) value = ReadBinary< Real >( file );
    for ( Real & value : z ) value = ReadBinary< Real >( file );

    for ( int i = 0; i < nodeCount; ++ i )
    {
        EXPECT_TRUE( std::isfinite( x[ i ] ) );
        EXPECT_TRUE( std::isfinite( y[ i ] ) );
        EXPECT_DOUBLE_EQ( z[ i ], 0.0 );
    }

    for ( int j = 0; j < 4; ++ j )
    {
        for ( int i = 0; i < 4; ++ i )
        {
            const int node = j * 4 + i;
            EXPECT_NEAR( x[ node ], static_cast< Real >( i ) / 3.0, 1.0e-12 );
            EXPECT_NEAR( y[ node ], static_cast< Real >( j ) / 3.0, 1.0e-12 );
        }
    }

    Real totalCellArea = 0.0;
    for ( int j = 0; j < 3; ++ j )
    {
        for ( int i = 0; i < 3; ++ i )
        {
            const int lowerLeft = j * 4 + i;
            const int lowerRight = lowerLeft + 1;
            const int upperLeft = lowerLeft + 4;
            const int upperRight = upperLeft + 1;
            const Real cellArea =
                0.5 * (
                    x[ lowerLeft ] * y[ lowerRight ] -
                    y[ lowerLeft ] * x[ lowerRight ] +
                    x[ lowerRight ] * y[ upperRight ] -
                    y[ lowerRight ] * x[ upperRight ] +
                    x[ upperRight ] * y[ upperLeft ] -
                    y[ upperRight ] * x[ upperLeft ] +
                    x[ upperLeft ] * y[ lowerLeft ] -
                    y[ upperLeft ] * x[ lowerLeft ] );

            EXPECT_GT( cellArea, 0.0 );
            EXPECT_NEAR( cellArea, 1.0 / 9.0, 1.0e-12 );
            totalCellArea += cellArea;
        }
    }

    EXPECT_NEAR( totalCellArea, 1.0, 1.0e-12 );
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

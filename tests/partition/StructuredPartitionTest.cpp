#include "StructuredPartition.h"

#include <gtest/gtest.h>

#include <cmath>
#include <numeric>
#include <vector>

using namespace ONEFLOW;

namespace
{

struct Grid3D
{
    int nx = 1;
    int ny = 1;
    int nz = 1;
    RealField x;
    RealField y;
    RealField z;

    void Build()
    {
        const int n = nx * ny * nz;
        x.resize( n );
        y.resize( n );
        z.resize( n );
        int id = 0;
        for ( int k = 0; k < nz; ++ k )
            for ( int j = 0; j < ny; ++ j )
                for ( int i = 0; i < nx; ++ i, ++ id )
                {
                    x[ id ] = static_cast< Real >( i );
                    y[ id ] = static_cast< Real >( j );
                    z[ id ] = static_cast< Real >( k );
                }
    }
};

std::vector< int > ZoneCounts( const IntField & zone, int nZones )
{
    std::vector< int > counts( static_cast< std::size_t >( nZones ), 0 );
    for ( int value : zone )
    {
        if ( value < 0 || value >= nZones ) continue;
        ++ counts[ value ];
    }
    return counts;
}

} // namespace

TEST( StructuredSplit, EffectiveLengthFollowsCellCountNotPhysicalSize )
{
    // Unit cube refined 40x40x160: the physical extent is 1 in every direction,
    // so only the median spacing can tell the axes apart.
    Grid3D grid;
    grid.nx = 40;
    grid.ny = 40;
    grid.nz = 160;
    grid.Build();
    for ( Real & value : grid.x ) value /= 40.0;
    for ( Real & value : grid.y ) value /= 40.0;
    for ( Real & value : grid.z ) value /= 160.0;

    const Real lx = EffectiveAxisLength( grid.x.data(), 40 * 40 * 160 );
    const Real ly = EffectiveAxisLength( grid.y.data(), 40 * 40 * 160 );
    const Real lz = EffectiveAxisLength( grid.z.data(), 40 * 40 * 160 );
    EXPECT_NEAR( lx, 39.0, 1.0e-9 );
    EXPECT_NEAR( ly, 39.0, 1.0e-9 );
    EXPECT_NEAR( lz, 159.0, 1.0e-9 );
    ( void ) grid;

    const StructuredSplit split = ChooseStructuredSplit( 32, lx, ly, lz );
    EXPECT_EQ( split.nx, 2 );
    EXPECT_EQ( split.ny, 2 );
    EXPECT_EQ( split.nz, 8 );
}

TEST( StructuredSplit, CubePrefersCubicBlocks )
{
    const StructuredSplit split = ChooseStructuredSplit( 8, 1.0, 1.0, 1.0 );
    EXPECT_EQ( split.nx, 2 );
    EXPECT_EQ( split.ny, 2 );
    EXPECT_EQ( split.nz, 2 );
}

TEST( StructuredSplit, AnisotropicBoxFollowsTheLongAxis )
{
    const StructuredSplit a = ChooseStructuredSplit( 16, 1.0, 1.0, 2.0 );
    EXPECT_EQ( a.nx * a.ny * a.nz, 16 );
    EXPECT_EQ( a.nz, 4 );

    const StructuredSplit b = ChooseStructuredSplit( 32, 1.0, 1.0, 4.0 );
    EXPECT_EQ( b.nx * b.ny * b.nz, 32 );
    EXPECT_EQ( b.nz, 8 );
}

TEST( StructuredSplit, LargeCountStaysCompact )
{
    const StructuredSplit split = ChooseStructuredSplit( 64, 1.0, 1.0, 1.0 );
    EXPECT_EQ( split.nx * split.ny * split.nz, 64 );
    EXPECT_EQ( split.nx, 4 );
    EXPECT_EQ( split.ny, 4 );
    EXPECT_EQ( split.nz, 4 );
}

TEST( StructuredSplit, SinglePartAndDegenerateExtents )
{
    const StructuredSplit one = ChooseStructuredSplit( 1, 1.0, 1.0, 1.0 );
    EXPECT_EQ( one.nx, 1 );
    EXPECT_EQ( one.ny, 1 );
    EXPECT_EQ( one.nz, 1 );

    // A zero thickness direction is never split.
    const StructuredSplit plate = ChooseStructuredSplit( 8, 1.0, 1.0, 0.0 );
    EXPECT_EQ( plate.nx * plate.ny * plate.nz, 8 );
    EXPECT_EQ( plate.nz, 1 );
}

TEST( StructuredZones, UniformCubeIsExactlyBalanced )
{
    Grid3D grid;
    grid.nx = 40;
    grid.ny = 40;
    grid.nz = 40;
    grid.Build();

    const StructuredSplit split = ChooseStructuredSplit( 8, 39.0, 39.0, 39.0 );
    IntField zone;
    ASSERT_TRUE( AssignStructuredZones(
        grid.nx * grid.ny * grid.nz, grid.x.data(), grid.y.data(),
        grid.z.data(), split, zone ) );

    const std::vector< int > counts = ZoneCounts( zone, 8 );
    for ( int count : counts ) EXPECT_EQ( count, 8000 );
}

TEST( StructuredZones, GradedMeshKeepsTheCellCountPerPart )
{
    // A non-uniform coordinate distribution must not unbalance the parts.
    Grid3D grid;
    grid.nx = 40;
    grid.ny = 1;
    grid.nz = 1;
    grid.Build();
    for ( int i = 0; i < grid.nx; ++ i )
    {
        grid.x[ i ] = static_cast< Real >( i ) * static_cast< Real >( i );
    }

    const StructuredSplit split = ChooseStructuredSplit( 8, 1521.0, 0.0, 0.0 );
    IntField zone;
    ASSERT_TRUE( AssignStructuredZones(
        grid.nx, grid.x.data(), grid.y.data(), grid.z.data(), split, zone ) );

    const std::vector< int > counts = ZoneCounts( zone, 8 );
    for ( int count : counts ) EXPECT_EQ( count, 5 );
}

TEST( StructuredZones, ZoneIdsStayInRangeAndCoverEveryPart )
{
    Grid3D grid;
    grid.nx = 30;
    grid.ny = 30;
    grid.nz = 30;
    grid.Build();

    const StructuredSplit split = ChooseStructuredSplit( 27, 29.0, 29.0, 29.0 );
    IntField zone;
    ASSERT_TRUE( AssignStructuredZones(
        grid.nx * grid.ny * grid.nz, grid.x.data(), grid.y.data(),
        grid.z.data(), split, zone ) );

    const std::vector< int > counts = ZoneCounts( zone, 27 );
    for ( int count : counts ) EXPECT_EQ( count, 1000 );
}

TEST( StructuredZones, RejectsInvalidInput )
{
    Grid3D grid;
    grid.nx = 4;
    grid.ny = 4;
    grid.nz = 4;
    grid.Build();

    const StructuredSplit split = ChooseStructuredSplit( 8, 3.0, 3.0, 3.0 );
    IntField zone;
    EXPECT_FALSE( AssignStructuredZones(
        0, grid.x.data(), grid.y.data(), grid.z.data(), split, zone ) );
    EXPECT_FALSE( AssignStructuredZones(
        64, nullptr, grid.y.data(), grid.z.data(), split, zone ) );

    StructuredSplit tooMany;
    tooMany.nx = 100;
    EXPECT_FALSE( AssignStructuredZones(
        64, grid.x.data(), grid.y.data(), grid.z.data(), tooMany, zone ) );
}

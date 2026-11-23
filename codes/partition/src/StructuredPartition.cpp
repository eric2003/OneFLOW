/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.
\*---------------------------------------------------------------------------*/

#include "StructuredPartition.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

BeginNameSpace( ONEFLOW )

namespace
{

bool HasExtent( Real value )
{
    return std::isfinite( value ) && value > 0.0;
}

// bin = rank in the coordinate order, scaled to nBins bins.  Ties are broken by
// the cell index so the result is independent of the sort implementation.
void AssignAxisBins(
    const Real * coordinate,
    int nCells,
    int nBins,
    IntField & bins )
{
    bins.assign( static_cast< std::size_t >( nCells ), 0 );
    if ( nBins <= 1 ) return;

    std::vector< int > order( static_cast< std::size_t >( nCells ) );
    std::iota( order.begin(), order.end(), 0 );
    std::stable_sort(
        order.begin(), order.end(),
        [ coordinate ]( int left, int right )
        {
            if ( coordinate[ left ] != coordinate[ right ] )
            {
                return coordinate[ left ] < coordinate[ right ];
            }
            return left < right;
        } );

    for ( int position = 0; position < nCells; ++ position )
    {
        const long long scaled =
            static_cast< long long >( position ) * nBins / nCells;
        bins[ order[ position ] ] = static_cast< int >( scaled );
    }
}

} // namespace

Real EffectiveAxisLength( const Real * coordinate, int nCells )
{
    if ( coordinate == nullptr || nCells < 2 ) return 0.0;

    Real lower = coordinate[ 0 ];
    Real upper = coordinate[ 0 ];
    for ( int cell = 1; cell < nCells; ++ cell )
    {
        lower = std::min( lower, coordinate[ cell ] );
        upper = std::max( upper, coordinate[ cell ] );
    }
    if ( ! std::isfinite( lower ) || ! std::isfinite( upper ) ) return 0.0;
    if ( upper <= lower ) return 0.0;

    std::vector< Real > sorted(
        coordinate, coordinate + static_cast< std::size_t >( nCells ) );
    std::sort( sorted.begin(), sorted.end() );

    // Cell centers that should share a plane differ by round-off, so a plain
    // "positive gap" test would report a spacing of ~1e-16 and make the
    // effective length explode.  Ignore gaps below a relative tolerance.
    const Real extent = upper - lower;
    const Real tolerance = extent * 1.0e-9;

    std::vector< Real > gaps;
    gaps.reserve( sorted.size() );
    for ( std::size_t i = 1; i < sorted.size(); ++ i )
    {
        const Real gap = sorted[ i ] - sorted[ i - 1 ];
        if ( gap > tolerance ) gaps.push_back( gap );
    }
    if ( gaps.empty() ) return 0.0;

    std::sort( gaps.begin(), gaps.end() );
    const Real median = gaps[ gaps.size() / 2 ];
    if ( median <= 0.0 ) return 0.0;
    return extent / median;
}

StructuredSplit ChooseStructuredSplit(
    int nPart, Real lx, Real ly, Real lz )
{
    StructuredSplit best;
    if ( nPart < 1 ) return best;

    const bool splitX = HasExtent( lx );
    const bool splitY = HasExtent( ly );
    const bool splitZ = HasExtent( lz );
    const Real inf = std::numeric_limits< Real >::infinity();

    best.nx = nPart;
    best.ny = 1;
    best.nz = 1;
    if ( ! splitX )
    {
        best.nx = 1;
        best.ny = nPart;
    }
    if ( ! splitX && ! splitY )
    {
        best.ny = 1;
        best.nz = nPart;
    }
    Real bestCost = inf;

    for ( int nx = 1; nx <= nPart; ++ nx )
    {
        if ( nPart % nx != 0 ) continue;
        if ( nx > 1 && ! splitX ) continue;
        const int rest = nPart / nx;
        for ( int ny = 1; ny <= rest; ++ ny )
        {
            if ( rest % ny != 0 ) continue;
            if ( ny > 1 && ! splitY ) continue;
            const int nz = rest / ny;
            if ( nz > 1 && ! splitZ ) continue;

            // Total internal interface area divided by the box volume:
            //   (nx-1)/lx + (ny-1)/ly + (nz-1)/lz
            Real cost = 0.0;
            if ( nx > 1 ) cost += static_cast< Real >( nx - 1 ) / lx;
            if ( ny > 1 ) cost += static_cast< Real >( ny - 1 ) / ly;
            if ( nz > 1 ) cost += static_cast< Real >( nz - 1 ) / lz;
            if ( cost < bestCost )
            {
                bestCost = cost;
                best.nx = nx;
                best.ny = ny;
                best.nz = nz;
            }
        }
    }
    return best;
}

bool AssignStructuredZones(
    int nCells,
    const Real * x,
    const Real * y,
    const Real * z,
    const StructuredSplit & split,
    IntField & zoneOfCell )
{
    if ( nCells <= 0 || x == nullptr || y == nullptr || z == nullptr )
        return false;
    if ( split.nx < 1 || split.ny < 1 || split.nz < 1 ) return false;
    if ( split.nx > nCells || split.ny > nCells || split.nz > nCells )
        return false;

    IntField binX;
    IntField binY;
    IntField binZ;
    AssignAxisBins( x, nCells, split.nx, binX );
    AssignAxisBins( y, nCells, split.ny, binY );
    AssignAxisBins( z, nCells, split.nz, binZ );

    zoneOfCell.resize( static_cast< std::size_t >( nCells ) );
    for ( int cell = 0; cell < nCells; ++ cell )
    {
        zoneOfCell[ cell ] = binX[ cell ]
            + split.nx * ( binY[ cell ] + split.ny * binZ[ cell ] );
    }
    return true;
}

EndNameSpace

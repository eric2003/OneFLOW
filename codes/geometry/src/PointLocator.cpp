/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "PointLocator.h"
#include <memory>
#include "GridHandles.h"
#include "Grid.h"
#include "NodeMesh.h"
#include "Constant.h"
#include "HXMath.h"
#include "Fatal.h"
#include <iostream>


BeginNameSpace( ONEFLOW )

class PointSpatialIndex
{
public:
    virtual ~PointSpatialIndex() = default;
    virtual void Insert( const RealField & coordinate, int pointId ) = 0;
    virtual void QueryRegion( const Real * pmin, const Real * pmax,
                              std::vector< int > & pointIds ) const = 0;
};

namespace
{
class AdtPointSpatialIndex final : public PointSpatialIndex
{
public:
    explicit AdtPointSpatialIndex( std::unique_ptr< AdtTree > tree )
        : tree( std::move( tree ) )
    {
    }

    void Insert( const RealField & coordinate, int pointId ) override
    {
        this->tree->AddNode(
            std::make_unique< AdtNode >( 3, coordinate.data(), pointId ) );
    }

    void QueryRegion( const Real * pmin, const Real * pmax,
                      std::vector< int > & pointIds ) const override
    {
        AdtTree::AdtNodeList nodes;
        this->tree->FindNodesInRegion( pmin, pmax, nodes );
        pointIds.reserve( pointIds.size() + nodes.size() );
        for ( const auto * node : nodes )
        {
            pointIds.push_back( node->GetData() );
        }
    }

private:
    std::unique_ptr< AdtTree > tree;
};
} // namespace

// Helper function to calculate squared distance (avoids expensive sqrt operation)
static inline Real CalcSquaredDistance(Real x1, Real y1, Real z1, Real x2, Real y2, Real z2) {
    Real dx = x1 - x2;
    Real dy = y1 - y2;
    Real dz = z1 - z2;
    return dx * dx + dy * dy + dz * dz;
}


PointLocator::PointLocator()
{
    this->id = 0;
    this->tolerance = 1.0e-6;
}

// [Refactored] unique_ptr automatically deletes the owned AdtTree.
PointLocator::~PointLocator() = default;


void PointLocator::Initialize( const RealField & pmin, const RealField & pmax, Real toleranceIn )
{
    this->tolerance = toleranceIn;
    std::unique_ptr< AdtTree > tree;
    ONEFLOW::CreateStandardADT( pmin, pmax, tree, this->tolerance );
    this->spatialIndex = std::make_unique< AdtPointSpatialIndex >( std::move( tree ) );
    this->ResetPointStore();
}

void PointLocator::Initialize( Grid * grid )
{
    std::unique_ptr< AdtTree > tree;
    ONEFLOW::CreateStandardADT( grid, tree, this->tolerance );
    this->spatialIndex = std::make_unique< AdtPointSpatialIndex >( std::move( tree ) );
    this->ResetPointStore();
}

void PointLocator::InitializeSpecial( Grid * grid, Real toleranceIn )
{
    // The ADT bounds must use the same tolerance as the narrow-phase query.
    // Reusing Initialize(grid) would build bounds with mindis / 4, then only
    // replace the query tolerance, which can exclude valid inserted points.
    grid->nodeMesh->CalcMinMaxBox();
    const RealField pmin = grid->nodeMesh->pmin;
    const RealField pmax = grid->nodeMesh->pmax;
    this->Initialize( pmin, pmax, toleranceIn );
}

void PointLocator::Initialize( Grids & grids )
{
    GridViews views = AsGridViews( grids );
    this->Initialize( views );
}

void PointLocator::Initialize( const GridViews & grids )
{
    std::unique_ptr< AdtTree > tree;
    ONEFLOW::CreateStandardADT( grids, tree, tolerance );
    this->spatialIndex = std::make_unique< AdtPointSpatialIndex >( std::move( tree ) );
    this->ResetPointStore();
}


void PointLocator::ResetPointStore()
{
    this->id = INVALID_INDEX;
    this->xCoor.clear();
    this->yCoor.clear();
    this->zCoor.clear();
}

void PointLocator::GetPoint( int id, Real & xm, Real & ym, Real & zm ) const
{
    xm = this->xCoor[ id ];
    ym = this->yCoor[ id ];
    zm = this->zCoor[ id ];
}

int PointLocator::AddPoint( Real xm, Real ym, Real zm )
{
    RealField coor( 3 );
    coor[ 0 ] = xm;
    coor[ 1 ] = ym;
    coor[ 2 ] = zm;
    return this->AddPoint( coor );
}

int PointLocator::FindPoint( Real xm, Real ym, Real zm )
{
    RealField coor( 3 );
    coor[ 0 ] = xm;
    coor[ 1 ] = ym;
    coor[ 2 ] = zm;
    return this->FindPoint( coor );
}

int PointLocator::FindPoint( const RealField & coordinate )
{
    std::vector< int > candidateIds;

    // 1. Broad Phase: Define the tolerance bounding box
    Real minWindow[3] = {
        coordinate[0] - this->tolerance,
        coordinate[1] - this->tolerance,
        coordinate[2] - this->tolerance
    };
    Real maxWindow[3] = {
        coordinate[0] + this->tolerance,
        coordinate[1] + this->tolerance,
        coordinate[2] + this->tolerance
    };

    // Find all candidate nodes within the bounding box
    this->spatialIndex->QueryRegion( minWindow, maxWindow, candidateIds );

    if ( candidateIds.empty() )
    {
        return INVALID_INDEX;
    }

    // 2. Narrow Phase: Find the closest point among candidates
    Real minSquaredDist = this->tolerance * this->tolerance;
    int bestId = INVALID_INDEX;

    for ( int candidateId : candidateIds )
    {
        const Real nx = this->xCoor[ candidateId ];
        const Real ny = this->yCoor[ candidateId ];
        const Real nz = this->zCoor[ candidateId ];

        const Real distSq = CalcSquaredDistance(
            coordinate[ 0 ], coordinate[ 1 ], coordinate[ 2 ],
            nx, ny, nz );

        if ( distSq < minSquaredDist )
        {
            minSquaredDist = distSq;
            bestId = candidateId;
        }
    }

    this->id = bestId;
    return this->id; // Returns INVALID_INDEX if no point is within exact tolerance
}

int PointLocator::AddPoint( const RealField & coor )
{
    // Reuse FindPoint logic to check for existing points within tolerance
    int existingId = this->FindPoint( coor );
    if ( existingId != INVALID_INDEX )
    {
        return existingId; // Point already exists, return its ID (Deduplication)
    }

    // Narrow Phase failed, add as a new unique point
    int newId = static_cast<int>( this->xCoor.size() );

    this->spatialIndex->Insert( coor, newId );

    this->xCoor.push_back( coor[0] );
    this->yCoor.push_back( coor[1] );
    this->zCoor.push_back( coor[2] );

    this->id = newId;
    return this->id;
}

void PointLocator::GetFaceCoorList( const IntField & nodeId, RealField &xList, RealField &yList, RealField &zList ) const
{
    for ( int i = 0; i < nodeId.size(); ++ i )
    {
        int ip = nodeId[ i ];
        xList.push_back( this->xCoor[ ip ] );
        yList.push_back( this->yCoor[ ip ] );
        zList.push_back( this->zCoor[ ip ] );
    }
}

void CreateStandardADT( const RealField & ptmin, const RealField & ptmax, std::unique_ptr<AdtTree>& adtTree, Real & tolerance )
{
    RealField pmin = ptmin;
    RealField pmax = ptmax;

    ONEFLOW::ShiftMinMaxBox( pmin, pmax, two * tolerance );

    adtTree = std::make_unique<AdtTree>( 3, pmin, pmax );
}

void CreateStandardADT( Grid * grid, std::unique_ptr<AdtTree>& adtTree, Real & tolerance )
{
    grid->nodeMesh->CalcMinMaxBox();
    RealField & ptmin = grid->nodeMesh->pmin;
    RealField & ptmax = grid->nodeMesh->pmax;

    RealField pmin = ptmin;
    RealField pmax = ptmax;

    Real mindis, maxdis;
    grid->GetMinMaxDistance( mindis, maxdis );

    tolerance = mindis / 4;

    ONEFLOW::ShiftMinMaxBox( pmin, pmax, two * tolerance );

    adtTree = std::make_unique<AdtTree>( 3, pmin, pmax );
}

void CreateStandardADT( Grids & grids, std::unique_ptr<AdtTree>& adtTree, Real & tolerance )
{
    GridViews views = AsGridViews( grids );
    CreateStandardADT( views, adtTree, tolerance );
}

void CreateStandardADT( const GridViews & grids, std::unique_ptr<AdtTree>& adtTree, Real & tolerance )
{
    RealField pmin( 3 ), pmax( 3 );
    ONEFLOW::GetBoundingBoxOfMultiZoneGrids( grids, pmin, pmax );

    tolerance = ONEFLOW::CalcGridTolerance( grids );
    if ( tolerance < 1.0e-10 )
    {
        tolerance = 1.0e-8;
    }

    ONEFLOW::ShiftMinMaxBox( pmin, pmax, two * tolerance );

    adtTree = std::make_unique<AdtTree>( 3, pmin, pmax );
}


void CreateStandardADTByTolerance( Grids & grids, std::unique_ptr<AdtTree>& adtTree, Real & tolerance )
{
    GridViews views = AsGridViews( grids );
    CreateStandardADTByTolerance( views, adtTree, tolerance );
}

void CreateStandardADTByTolerance( const GridViews & grids, std::unique_ptr<AdtTree>& adtTree, Real & tolerance )
{
    RealField pmin( 3 ), pmax( 3 );
    ONEFLOW::GetBoundingBoxOfMultiZoneGrids( grids, pmin, pmax );

    ONEFLOW::ShiftMinMaxBox( pmin, pmax, two * tolerance );

    adtTree = std::make_unique<AdtTree>( 3, pmin, pmax );
}

void ShiftMinMaxBox( RealField & pmin, RealField & pmax, Real tolerance )
{
    pmin[ 0 ] -= tolerance;
    pmin[ 1 ] -= tolerance;
    pmin[ 2 ] -= tolerance;

    pmax[ 0 ] += tolerance;
    pmax[ 1 ] += tolerance;
    pmax[ 2 ] += tolerance;
}

void GetGridsMinMaxDistance( Grids & grids, Real & mindis, Real & maxdis )
{
    GridViews views = AsGridViews( grids );
    GetGridsMinMaxDistance( views, mindis, maxdis );
}

void GetGridsMinMaxDistance( const GridViews & grids, Real & mindis, Real & maxdis )
{
    mindis =   LARGE;
    maxdis = - LARGE;

    int numberOfZones = static_cast< int >( grids.size() );

    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        Real dismin, dismax;
        GridAt( grids, iZone )->GetMinMaxDistance( dismin, dismax );

        mindis = ONEFLOW::MIN( mindis, dismin );
        maxdis = ONEFLOW::MAX( maxdis, dismax );
    }
}

Real CalcGridTolerance( Grids & grids )
{
    GridViews views = AsGridViews( grids );
    return CalcGridTolerance( views );
}

Real CalcGridTolerance( const GridViews & grids )
{
    Real mindis =   LARGE;
    Real maxdis = - LARGE;

    ONEFLOW::GetGridsMinMaxDistance( grids, mindis, maxdis );

    Real tolerance = mindis / 10.0;

    return tolerance;
}

void GetBoundingBoxOfMultiZoneGrids( Grids & grids, RealField & pmin, RealField & pmax )
{
    GridViews views = AsGridViews( grids );
    GetBoundingBoxOfMultiZoneGrids( views, pmin, pmax );
}

void GetBoundingBoxOfMultiZoneGrids( const GridViews & grids, RealField & pmin, RealField & pmax )
{
    pmin[ 0 ] = LARGE;
    pmin[ 1 ] = LARGE;
    pmin[ 2 ] = LARGE;

    pmax[ 0 ] = - LARGE;
    pmax[ 1 ] = - LARGE;
    pmax[ 2 ] = - LARGE;

    int numberOfZones = static_cast< int >( grids.size() );

    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        GridAt( grids, iZone )->nodeMesh->CalcMinMaxBox();
        RealField & localPmin = GridAt( grids, iZone )->nodeMesh->pmin;
        RealField & localPmax = GridAt( grids, iZone )->nodeMesh->pmax;

        pmin[ 0 ] = ONEFLOW::MIN( pmin[ 0 ], localPmin[ 0 ] );
        pmin[ 1 ] = ONEFLOW::MIN( pmin[ 1 ], localPmin[ 1 ] );
        pmin[ 2 ] = ONEFLOW::MIN( pmin[ 2 ], localPmin[ 2 ] );

        pmax[ 0 ] = ONEFLOW::MAX( pmax[ 0 ], localPmax[ 0 ] );
        pmax[ 1 ] = ONEFLOW::MAX( pmax[ 1 ], localPmax[ 1 ] );
        pmax[ 2 ] = ONEFLOW::MAX( pmax[ 2 ], localPmax[ 2 ] );
    }
}

EndNameSpace

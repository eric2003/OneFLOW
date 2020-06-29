/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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

#include "PointSearch.h"
#include "Grid.h"
#include "NodeMesh.h"
#include "HXMath.h"
#include "Stop.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

PointSearch::PointSearch()
{
    this->coorTree = 0;
}

PointSearch::~PointSearch()
{
    delete this->coorTree;
}

void PointSearch::Initialize( RealField & pmin, RealField & pmax, Real toleranceIn )
{
    this->tolerance = toleranceIn;
    ONEFLOW::CreateStandardADT( pmin, pmax, this->coorTree, this->tolerance );
}

void PointSearch::Initialize( Grid * grid )
{
    ONEFLOW::CreateStandardADT( grid, this->coorTree, this->tolerance );
}

void PointSearch::InitializeSpecial( Grid * grid, Real toleranceIn )
{
    this->Initialize( grid );
    this->tolerance = toleranceIn;
}

void PointSearch::Initialize( Grids & grids )
{
    ONEFLOW::CreateStandardADT( grids, this->coorTree, tolerance );
}

int PointSearch::AddPoint( RealField & coor )
{
    Real minWindow[ 3 ];
    Real maxWindow[ 3 ];

    minWindow[ 0 ] = coor[ 0 ] - tolerance;
    minWindow[ 1 ] = coor[ 1 ] - tolerance;
    minWindow[ 2 ] = coor[ 2 ] - tolerance;

    maxWindow[ 0 ] = coor[ 0 ] + tolerance;
    maxWindow[ 1 ] = coor[ 1 ] + tolerance;
    maxWindow[ 2 ] = coor[ 2 ] + tolerance;

    AdtTree::AdtNodeList nodeList;
    this->coorTree->FindNodesInRegion( minWindow, maxWindow, nodeList );

    if ( nodeList.size() == 0 )
    {
        int count = this->xCoor.size();
        AdtNode * node = new AdtNode( 3, & coor[ 0 ], count );
        this->coorTree->AddNode( node );
        this->id = this->xCoor.size();

        xCoor.push_back( coor[ 0 ] );
        yCoor.push_back( coor[ 1 ] );
        zCoor.push_back( coor[ 2 ] );
        
        return this->id;
    }
    else
    {
        if ( nodeList.size() > 1 )
        {
            cout << "FATAL ERROR : nodeList.size() = " << nodeList.size()<<endl;
            Stop("");
        }
        AdtNode * node =  nodeList[ 0 ];

        this->id = node->GetData();

        return this->id;
    }
}

void PointSearch::GetPoint( int id, Real & xm, Real & ym, Real & zm )
{
    xm = this->xCoor[ id ];
    ym = this->yCoor[ id ];
    zm = this->zCoor[ id ];
}

int PointSearch::AddPoint( Real xm, Real ym, Real zm )
{
    RealField coor( 3 );
    coor[ 0 ] = xm;
    coor[ 1 ] = ym;
    coor[ 2 ] = zm;
    return this->AddPoint( coor );
}

int PointSearch::FindPoint( Real xm, Real ym, Real zm )
{
    RealField coor( 3 );
    coor[ 0 ] = xm;
    coor[ 1 ] = ym;
    coor[ 2 ] = zm;
    return this->FindPoint( coor );
}

int PointSearch::FindPoint( RealField & coordinate )
{
    AdtTree::AdtNodeList nodeList;

    Real minWindow[ 3 ];
    Real maxWindow[ 3 ];

    minWindow[ 0 ] = coordinate[ 0 ] - this->tolerance;
    minWindow[ 1 ] = coordinate[ 1 ] - this->tolerance;
    minWindow[ 2 ] = coordinate[ 2 ] - this->tolerance;

    maxWindow[ 0 ] = coordinate[ 0 ] + this->tolerance;
    maxWindow[ 1 ] = coordinate[ 1 ] + this->tolerance;
    maxWindow[ 2 ] = coordinate[ 2 ] + this->tolerance;

    nodeList.resize( 0 );
    this->coorTree->FindNodesInRegion( minWindow, maxWindow, nodeList );

    if ( nodeList.size() == 0 )
    {
        return INVALID_INDEX;
    }
    else
    {
        if ( nodeList.size() > 1 )
        {
            int numberOfSize = nodeList.size();
            cout << " impossible nodeList.size() = " << nodeList.size() << endl;
            
            int kkk = 1;
            Stop( "" );
        }
        AdtNode * node = nodeList[ 0 ];
        this->id = node->GetData();
        return this->id;
    }
}

void PointSearch::GetFaceCoorList( IntField & nodeId, RealField &xList, RealField &yList, RealField &zList )
{
    for ( int i = 0; i < nodeId.size(); ++ i )
    {
        int ip = nodeId[ i ];
        xList.push_back( this->xCoor[ ip ] );
        yList.push_back( this->yCoor[ ip ] );
        zList.push_back( this->zCoor[ ip ] );
    }
}

void CreateStandardADT( RealField & ptmin, RealField & ptmax, AdtTree *& adtTree, Real & tolerance )
{
    RealField pmin = ptmin;
    RealField pmax = ptmax;

    ONEFLOW::ShiftMinMaxBox( pmin, pmax, two * tolerance );

    adtTree = new AdtTree( 3, pmin, pmax );
}

void CreateStandardADT( Grid * grid, AdtTree *& adtTree, Real & tolerance )
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

    adtTree = new AdtTree( 3, pmin, pmax );
}

void CreateStandardADT( Grids & grids, AdtTree *& adtTree, Real & tolerance )
{
    RealField pmin( 3 ), pmax( 3 );
    ONEFLOW::GetBoundingBoxOfMultiZoneGrids( grids, pmin, pmax );

    tolerance = ONEFLOW::CalcGridTolerance( grids );
    if ( tolerance < 1.0e-10 )
    {
        tolerance = 1.0e-8;
    }

    ONEFLOW::ShiftMinMaxBox( pmin, pmax, two * tolerance );

    adtTree = new AdtTree( 3, pmin, pmax );
}


void CreateStandardADTByTolerance( Grids & grids, AdtTree *& adtTree, Real & tolerance )
{
    RealField pmin( 3 ), pmax( 3 );
    ONEFLOW::GetBoundingBoxOfMultiZoneGrids( grids, pmin, pmax );

    ONEFLOW::ShiftMinMaxBox( pmin, pmax, two * tolerance );

    adtTree = new AdtTree( 3, pmin, pmax );
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
    mindis =   LARGE;
    maxdis = - LARGE;

    int numberOfZones = grids.size();

    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        Real dismin, dismax;
        grids[ iZone ]->GetMinMaxDistance( dismin, dismax );

        mindis = ONEFLOW::MIN( mindis, dismin );
        maxdis = ONEFLOW::MAX( maxdis, dismax );
    }
}

Real CalcGridTolerance( Grids & grids )
{
    Real mindis =   LARGE;
    Real maxdis = - LARGE;

    ONEFLOW::GetGridsMinMaxDistance( grids, mindis, maxdis );

    Real tolerance = mindis / 10.0;

    return tolerance;
}

void GetBoundingBoxOfMultiZoneGrids( Grids & grids, RealField & pmin, RealField & pmax )
{
    pmin[ 0 ] = LARGE;
    pmin[ 1 ] = LARGE;
    pmin[ 2 ] = LARGE;

    pmax[ 0 ] = - LARGE;
    pmax[ 1 ] = - LARGE;
    pmax[ 2 ] = - LARGE;

    int numberOfZones = grids.size();

    for ( int iZone = 0; iZone < numberOfZones; ++ iZone )
    {
        grids[ iZone ]->nodeMesh->CalcMinMaxBox();
        RealField & localPmin = grids[ iZone ]->nodeMesh->pmin;
        RealField & localPmax = grids[ iZone ]->nodeMesh->pmax;

        pmin[ 0 ] = ONEFLOW::MIN( pmin[ 0 ], localPmin[ 0 ] );
        pmin[ 1 ] = ONEFLOW::MIN( pmin[ 1 ], localPmin[ 1 ] );
        pmin[ 2 ] = ONEFLOW::MIN( pmin[ 2 ], localPmin[ 2 ] );

        pmax[ 0 ] = ONEFLOW::MAX( pmax[ 0 ], localPmax[ 0 ] );
        pmax[ 1 ] = ONEFLOW::MAX( pmax[ 1 ], localPmax[ 1 ] );
        pmax[ 2 ] = ONEFLOW::MAX( pmax[ 2 ], localPmax[ 2 ] );
    }
}

EndNameSpace
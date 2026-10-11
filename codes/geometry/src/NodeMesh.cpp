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

#include "NodeMesh.h"
#include "HXMath.h"
#include <utility>
#include <stdexcept>

BeginNameSpace( ONEFLOW )

namespace
{
void ValidateCoordinateArrays( const NodeMesh & nodeMesh )
{
    if ( nodeMesh.xN.size() != nodeMesh.yN.size() ||
         nodeMesh.xN.size() != nodeMesh.zN.size() )
    {
        throw std::logic_error( "NodeMesh coordinate arrays have inconsistent sizes" );
    }
}
} // namespace

NodeMesh::NodeMesh()
{
    boxFlag = false;
    pmin.resize( 3 );
    pmax.resize( 3 );
}

NodeMesh::Point3 NodeMesh::GetPoint( HXSize_t nodeIndex ) const
{
    ValidateCoordinateArrays( *this );
    if ( nodeIndex >= this->xN.size() )
    {
        throw std::out_of_range( "NodeMesh node index is out of range" );
    }
    return { this->xN[ nodeIndex ], this->yN[ nodeIndex ], this->zN[ nodeIndex ] };
}

void NodeMesh::SetPoint( HXSize_t nodeIndex, Real xp, Real yp, Real zp )
{
    ValidateCoordinateArrays( *this );
    if ( nodeIndex >= this->xN.size() )
    {
        throw std::out_of_range( "NodeMesh node index is out of range" );
    }
    this->xN[ nodeIndex ] = xp;
    this->yN[ nodeIndex ] = yp;
    this->zN[ nodeIndex ] = zp;
    this->boxFlag = false;
}

void NodeMesh::CreateNodes( int numberOfNodes )
{
    if ( numberOfNodes < 0 )
    {
        throw std::invalid_argument( "NodeMesh node count cannot be negative" );
    }

    xN.resize( static_cast< HXSize_t >( numberOfNodes ) );
    yN.resize( static_cast< HXSize_t >( numberOfNodes ) );
    zN.resize( static_cast< HXSize_t >( numberOfNodes ) );
    this->boxFlag = false;
}

void NodeMesh::SetCoordinates( RealField x, RealField y, RealField z )
{
    if ( x.size() != y.size() || x.size() != z.size() )
    {
        throw std::invalid_argument(
            "NodeMesh coordinate arrays must have equal sizes" );
    }

    this->xN = std::move( x );
    this->yN = std::move( y );
    this->zN = std::move( z );
    this->boxFlag = false;
}

void NodeMesh::CalcMinMaxBox()
{
    // Coordinate arrays are public, so callers can mutate them without invalidating
    // the cached bounds. Recompute to keep the bounding box correct at this boundary.
    ValidateCoordinateArrays( *this );

    const HXSize_t numberOfNodes = this->GetNumberOfNodes();
    if ( numberOfNodes == 0 )
    {
        for ( int iDim = 0; iDim < 3; ++ iDim )
        {
            this->pmin[ iDim ] = 0;
            this->pmax[ iDim ] = 0;
        }
        this->boxFlag = true;
        return;
    }

    Real xmin =   LARGE;
    Real ymin =   LARGE;
    Real zmin =   LARGE;
    Real xmax = - LARGE;
    Real ymax = - LARGE;
    Real zmax = - LARGE;

    for ( HXSize_t iNode = 0; iNode < numberOfNodes; ++ iNode )
    {
        xmin = ONEFLOW::MIN( xmin, xN[ iNode ] );
        ymin = ONEFLOW::MIN( ymin, yN[ iNode ] );
        zmin = ONEFLOW::MIN( zmin, zN[ iNode ] );

        xmax = ONEFLOW::MAX( xmax, xN[ iNode ] );
        ymax = ONEFLOW::MAX( ymax, yN[ iNode ] );
        zmax = ONEFLOW::MAX( zmax, zN[ iNode ] );
    }

    pmin[ 0 ] = xmin;
    pmin[ 1 ] = ymin;
    pmin[ 2 ] = zmin;

    pmax[ 0 ] = xmax;
    pmax[ 1 ] = ymax;
    pmax[ 2 ] = zmax;
    this->boxFlag = true;
}

void NodeMesh::AddPoint( Real xp, Real yp, Real zp )
{
    ValidateCoordinateArrays( *this );
    this->xN.push_back( xp );
    this->yN.push_back( yp );
    this->zN.push_back( zp );
    this->boxFlag = false;
}

void NodeMesh::ScaleAndTranslate( Real scale, Real tx, Real ty, Real tz )
{
    ValidateCoordinateArrays( *this );
    const size_t nNodes = this->GetNumberOfNodes();
    for ( size_t iNode = 0; iNode < nNodes; ++ iNode )
    {
        this->xN[ iNode ] *= scale;
        this->yN[ iNode ] *= scale;
        this->zN[ iNode ] *= scale;

        this->xN[ iNode ] += tx;
        this->yN[ iNode ] += ty;
        this->zN[ iNode ] += tz;
    }
    this->boxFlag = false;
}

void NodeMesh::TurnZAxisToYAxis()
{
    ValidateCoordinateArrays( *this );
    const size_t nNodes = this->GetNumberOfNodes();
    for ( size_t iNode = 0; iNode < nNodes; ++ iNode )
    {
        const Real y = this->yN[ iNode ];
        this->yN[ iNode ] = this->zN[ iNode ];
        this->zN[ iNode ] = - y;
    }
    this->boxFlag = false;
}


EndNameSpace

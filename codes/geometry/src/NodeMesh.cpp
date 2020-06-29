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

#include "NodeMesh.h"
#include "HXMath.h"

BeginNameSpace( ONEFLOW )

NodeMesh::NodeMesh()
{
    boxFlag = false;
    pmin.resize( 3 );
    pmax.resize( 3 );
}

NodeMesh::NodeMesh( const NodeMesh & rhs )
{
    this->boxFlag = rhs.boxFlag;
    this->pmin = rhs.pmin;
    this->pmax = rhs.pmax;
    this->xN = rhs.xN;
    this->yN = rhs.yN;
    this->zN = rhs.zN;
}

NodeMesh::~NodeMesh()
{
    ;
}

void NodeMesh::CreateNodes( int numberOfNodes )
{
    xN.resize( numberOfNodes );
    yN.resize( numberOfNodes );
    zN.resize( numberOfNodes );
}

void NodeMesh::CalcMinMaxBox()
{
    if ( this->boxFlag ) return;

    int numberOfNodes = this->GetNumberOfNodes();

    Real xmin =   LARGE;
    Real ymin =   LARGE;
    Real zmin =   LARGE;
    Real xmax = - LARGE;
    Real ymax = - LARGE;
    Real zmax = - LARGE;

    for ( int iNode = 0; iNode < numberOfNodes; ++ iNode )
    {
        Real tmp = ABS( zN[ iNode ] );
        if ( tmp > 1.0e10 )
        {
            int kkk = 1;
        }
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
}

void NodeMesh::AddPoint( Real xp, Real yp, Real zp )
{
    this->xN.push_back( xp );
    this->yN.push_back( yp );
    this->zN.push_back( zp );
}


EndNameSpace
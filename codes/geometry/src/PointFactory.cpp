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

#include "PointFactory.h"
#include "Grid.h"
#include "NodeMesh.h"
#include "HXMath.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

PointBasic::PointBasic()
{
}

PointBasic::~PointBasic()
{
}

bool PointBasic::FindPoint( PointBasic::PointType & point, PointBasic::PointSet::iterator & iter )
{
    iter = this->pointSet.find( point );
    if ( iter == this->pointSet.end() )
    {
        return false;
    }
    return true;
}

int PointBasic::FindPoint( Real xm, Real ym, Real zm )
{
    PointBasic::PointType pt( xm, ym, zm );
    return this->FindPointId( pt );
}

int PointBasic::FindPointId( PointBasic::PointType & point )
{
    PointBasic::PointSet::iterator iter = this->pointSet.find( point );
    if ( iter == this->pointSet.end() )
    {
        return -1;
    }
    return iter->id;
}

int PointBasic::AddPoint( Real xm, Real ym, Real zm )
{
    PointBasic::PointType pt( xm, ym, zm );
    return this->AddPoint( pt );
}

int PointBasic::DeletePoint( Real xm, Real ym, Real zm )
{
    PointBasic::PointType pt( xm, ym, zm );
    return this->DeletePoint( pt );
}

int PointBasic::DeletePoint( PointBasic::PointType & point )
{
    PointSet::iterator iter;
    int pid = -1;
    if ( this->FindPoint( point, iter ) )
    {
        pid = iter->id;
        pointSet.erase( point );
        int kkk = 1;
    }

    return pid;
}

int PointBasic::AddPoint( PointBasic::PointType & point )
{
    PointSet::iterator iter;
    int pid;
    if ( this->FindPoint( point, iter ) )
    {
        pid = iter->id;
    }
    else
    {
        int index = this->pointList.size();
        point.id = index;
        this->AddPointDirectly( point );
        pid = index;
    }

    return pid;
}

void PointBasic::AddPointDirectly( PointBasic::PointType & point, int bcType )
{
    this->pointSet.insert( point );
    this->pointList.push_back( point );
}

void PointBasic::GetFaceCoorList( IntField & nodeId, RealField &xList, RealField &yList, RealField &zList )
{
    for ( int i = 0; i < nodeId.size(); ++ i )
    {
        int ip = nodeId[ i ];
        PointBasic::PointType & pt = pointList[ ip ];
        xList.push_back( pt.x );
        yList.push_back( pt.y );
        zList.push_back( pt.z );
    }
}


PointFactory::PointFactory()
{
}

PointFactory::~PointFactory()
{
}

void PointFactory::InitC2g()
{
    int nPoint = this->pointList.size();
    this->c2g.resize( nPoint );
    for ( int iNode = 0; iNode < nPoint; ++ iNode )
    {
        this->c2g[ iNode ] = iNode;
    }
}

EndNameSpace
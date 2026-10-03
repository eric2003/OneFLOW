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

#include "PointMachine.h"


BeginNameSpace( ONEFLOW )

PointMachine point_Machine;

PointMachine::PointMachine()
{
    ;
}

PointMachine::~PointMachine() = default;

void PointMachine::Reset()
{
    ptList.clear();
}

void PointMachine::AddPoint( Real x, Real y, Real z, int id )
{
    auto pt = std::make_unique< PointType >( x, y, z, id );
    this->ptList.push_back( std::move( pt ) );

}

PointType & PointMachine::GetPoint( int id )
{
    const int index = id - 1;
    return * this->ptList[ index ];
}

const PointType & PointMachine::GetPoint( int id ) const
{
    const int index = id - 1;
    return * this->ptList.at( index );
}

EndNameSpace

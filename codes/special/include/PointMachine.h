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


#pragma once
#include "HXDefine.h"
#include "Point.h"
#include <memory>

BeginNameSpace( ONEFLOW )

using PointType = Point< Real >;

class PointMachine
{
public:
    PointMachine();
    ~PointMachine();
public:
    void Reset();
    void AddPoint( Real x, Real y, Real z, int id = 0 );
    PointType & GetPoint( int id );
    const PointType & GetPoint( int id ) const;
    int GetNPoint() const { return static_cast< int >( this->ptList.size() ); };
private:
    // PointMachine owns the layout points; callers receive non-owning views.
    HXVector< std::unique_ptr< PointType > > ptList;
};

extern PointMachine point_Machine;
EndNameSpace

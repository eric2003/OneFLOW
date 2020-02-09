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

#include "PointMachine.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

PointMachine point_Machine;

PointMachine::PointMachine()
{
    ;
}

PointMachine::~PointMachine()
{
    for ( int i = 0; i < ptList.size(); ++ i )
    {
        delete ptList[ i ];
    }
}

void PointMachine::AddPoint( Real x, Real y, Real z, int id )
{
    PointType * pt = new PointType( x, y, z, id );
    this->ptList.push_back( pt );
    int idd = ptBasic.AddPoint( x, y, z );
    //int idd1 = ptBasic.DeletePoint( x, y, z );
    int kkk = 1;
}

PointType * PointMachine::GetPoint( int id )
{
    int ida = id - 1;
    PointType * pt = this->ptList[ ida ];
    return pt;
}

EndNameSpace
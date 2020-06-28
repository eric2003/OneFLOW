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

#include "CalcCoor.h"

BeginNameSpace( ONEFLOW )

CalcCoor::CalcCoor()
{
}

CalcCoor::~CalcCoor()
{
}

void CalcCoor::SetCoor( int i, int j, int k )
{
    this->i = i;
    this->j = j;
    this->k = k;
}

bool CalcCoor::operator < ( const CalcCoor & rhs ) const
{
    if ( this->i == rhs.i )
    {
        if ( this->j == rhs.j )
        {
            return ( this->k < rhs.k );
        }
        else
        {
            return ( this->j < rhs.j );
        }
    }
    else
    {
        return ( this->i < rhs.i );
    }
}

EndNameSpace
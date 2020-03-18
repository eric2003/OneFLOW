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
#pragma once
#include "HXVector.h"

BeginNameSpace( ONEFLOW )

template < typename T >
class Marray
{
public:
    Marray()
    {
    }

    Marray( UInt nEqu, int numberOfCells )
    {
        data.resize( nEqu );
        for ( UInt iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            data[ iEqu ].resize( numberOfCells );
        }
    }
    ~Marray()
    {
    }
protected:
    HXVector< HXVector< T > > data;
public:
    UInt GetNEqu() { return data.size(); }

    HXVector< T > & operator[]( int iEqu )
    {
        return data[ iEqu ];
    }

    Marray< T > & operator = ( const T & value )
    {
        UInt nEqu = this->GetNEqu();
        for ( UInt iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            data[ iEqu ] = value;
        }
        return * this;
    }
};


EndNameSpace

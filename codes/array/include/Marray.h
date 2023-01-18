/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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
#include "HXType.h"
#include "HXVector.h"

BeginNameSpace( ONEFLOW )

template < typename T >
class Marray
{
public:
    Marray()
    {
    }

    Marray( HXSize_t nEqu, int numberOfCells )
    {
        data.resize( nEqu );
        for ( HXSize_t iEqu = 0; iEqu < nEqu; ++ iEqu )
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
    HXSize_t GetNEqu() { return data.size(); }

    HXVector< T > & operator[]( int iEqu )
    {
        return data[ iEqu ];
    }

    HXVector< T > & AsOneD()
    {
        return data[ 0 ];
    }

    Marray< T > & operator = ( const T & value )
    {
        HXSize_t nEqu = this->GetNEqu();
        for ( HXSize_t iEqu = 0; iEqu < nEqu; ++ iEqu )
        {
            data[ iEqu ] = value;
        }
        return * this;
    }
};


EndNameSpace

/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
	Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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

#include "Mid.h"

BeginNameSpace( ONEFLOW )

Mid::Mid()
{
    this->size = 0;
    this->id   = 0;
}

Mid::Mid( int size, int id )
{
    this->size = size;
    this->id   = id;
	this->data.resize( size );
}

Mid::Mid( const Mid & rhs )
{
    this->size = rhs.size;
    this->id   = rhs.id;
	this->data.resize( size );

    for ( int i = 0; i < size; ++ i )
    {
        this->data[ i ] = rhs.data[ i ];
    }
}

Mid & Mid::operator = ( const Mid & rhs )
{
    if ( this == & rhs ) return * this;

    if ( this->size != rhs.size )
    {
        this->size = rhs.size;
		this->data.resize( this->size );
    }

    this->id   = rhs.id;
    for ( int i = 0; i < size; ++ i )
    {
        this->data[ i ] = rhs.data[ i ];
    }

    return * this;
}

Mid::~Mid()
{
}

bool Mid::operator < ( const Mid & rhs ) const
{
    if ( this->size != rhs.size ) return this->size < rhs.size;
    for ( int i = 0; i < this->size; ++ i )
    {
        if ( this->data[ i ] != rhs.data[ i ] )
        {
            return this->data[ i ] < rhs.data[ i ];
        }
    }
    return false;
}


EndNameSpace
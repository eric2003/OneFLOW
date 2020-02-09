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
#include "Configure.h"

BeginNameSpace( ONEFLOW )

class Range
{
public:
    Range()
    {
        first  = 0;
        last   = 0;
        stride = 1;
    }

    Range( const Range& r )
    {
        first  = r.first;
        last   = r.last;
        stride = r.stride;
    }

    explicit Range( int slicePosition )
    {
        first  = slicePosition;
        last   = slicePosition;
        stride = 1;
    }

    Range( int firstIn, int lastIn, int strideIn = 1 )
        : first( firstIn ), last( lastIn ), stride( strideIn )
    { 
    }

    int First() const
    { 
        return first; 
    }

    int Last() const
    {
        return last;
    }

    unsigned Length() const
    {
        return ( last - first ) / stride + 1;
    }

    int Stride() const
    {
        return stride;
    }

    void SetRange( int first, int last, int stride = 1 )
    {
        this->first  = first;
        this->last   = last;
        this->stride = stride;
    }

    Range operator - ( int shift ) const
    { 
        return Range( first - shift, last - shift, stride ); 
    }

    Range operator + ( int shift ) const
    { 
        return Range( first + shift, last + shift, stride ); 
    }

    int operator [] ( unsigned i ) const
    {
        return first + i * stride;
    }

    int operator () ( unsigned i ) const
    {
        return first + i * stride;
    }
protected:
    int first, last, stride;
};

EndNameSpace
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
#include "NamespaceMacros.h"
#include "Range.h"
#include <memory>

BeginNameSpace( ONEFLOW )

template < typename T >
class ArrayPointer;

template < typename T >
class ArrayPointer
{
public:
    ArrayPointer()
        : shareMem( nullptr )
        , data1( nullptr )
        , data2( nullptr )
        , data3( nullptr )
        , data4( nullptr )
        , datap1( nullptr )
        , datap2( nullptr )
        , datap3( nullptr )
        , datap4( nullptr )
    {
    }

    // unique_ptr members release owned blocks; shareMem views are not owned.
    ~ArrayPointer() = default;

public:
    T * shareMem;
    // Non-owning views used by Multiarray indexing (may point into shareMem).
    T *    data1;
    T **   data2;
    T ***  data3;
    T **** data4;

    T *    datap1;
    T **   datap2;
    T ***  datap3;
    T **** datap4;

private:
    // Owned storage for index tables and element buffers.
    std::unique_ptr<T[]>    owned1;
    std::unique_ptr<T*[]>   owned2;
    std::unique_ptr<T**[]>  owned3;
    std::unique_ptr<T***[]> owned4;

public:

    void Allocate( Range r0 )
    {
        this->AllocateArray( r0 );
    }

    void Allocate( Range r0, Range r1 )
    {
        this->AllocateArray( r0, r1 );
    }

    void Allocate( Range r0, Range r1, Range r2 )
    {
        this->AllocateArray( r0, r1, r2 );
    }

    void Allocate( Range r0, Range r1, Range r2, Range r3 )
    {
        this->AllocateArray( r0, r1, r2, r3 );
    }

    void Allocate( T * objects, Range r0 )
    {
        this->shareMem = objects;

        this->AllocateArray( r0 );
    }

    void Allocate( T * objects, Range r0, Range r1 )
    {
        this->shareMem = objects;

        this->AllocateArray( r0, r1 );
    }

    void Allocate( T * objects, Range r0, Range r1, Range r2 )
    {
        this->shareMem = objects;

        this->AllocateArray( r0, r1, r2 );
    }

    void Allocate( T * objects, Range r0, Range r1, Range r2, Range r3 )
    {
        this->shareMem = objects;

        this->AllocateArray( r0, r1, r2, r3 );
    }

    void AllocateArray( Range r0, T * objects )
    {
        int bound0 = r0.Length();

        if ( ! shareMem )
        {
            this->owned1 = std::make_unique<T[]>( static_cast<std::size_t>( bound0 ) );
            this->data1 = this->owned1.get();
        }
        else
        {
            this->owned1.reset();
            this->data1 = shareMem;
        }

        int st1 = r0.First();

        this->datap1 = this->data1 - st1;
    }

    void AllocateArray( Range r0 )
    {
        this->AllocateArray( r0, static_cast<T *>( nullptr ) );
    }

    void AllocateArray( Range r0, Range r1 )
    {
        int bound0 = r0.Length();
        int bound1 = r1.Length();

        int bound01 = bound0 * bound1;
        this->owned2 = std::make_unique<T*[]>( static_cast<std::size_t>( bound1 ) );
        this->data2 = this->owned2.get();

        if ( ! shareMem )
        {
            this->owned1 = std::make_unique<T[]>( static_cast<std::size_t>( bound01 ) );
            this->data1 = this->owned1.get();
        }
        else
        {
            this->owned1.reset();
            this->data1 = shareMem;
        }

        int st1 = r0.First();
        int st2 = r1.First();

        this->datap2 = this->data2 - st2;

        for ( int i1 = 1; i1 < bound1; ++ i1 )
        {
            int j0 = bound0 * i1;  // = bound0*(i1 + j1) where j1 = 0
            this->datap2[ i1 + st2 ] = & this->data1[ j0 ] - st1;
        }
    }

    void AllocateArray( Range r0, Range r1, Range r2 )
    {
        int bound0 = r0.Length();
        int bound1 = r1.Length();
        int bound2 = r2.Length();

        int bound12  = bound1 * bound2;
        int bound012 = bound0 * bound12;
        this->owned3 = std::make_unique<T**[]>( static_cast<std::size_t>( bound2 ) );
        this->owned2 = std::make_unique<T*[]>( static_cast<std::size_t>( bound12 ) );
        this->data3 = this->owned3.get();
        this->data2 = this->owned2.get();

        if ( ! shareMem )
        {
            this->owned1 = std::make_unique<T[]>( static_cast<std::size_t>( bound012 ) );
            this->data1 = this->owned1.get();
        }
        else
        {
            this->owned1.reset();
            this->data1 = shareMem;
        }

        int st1 = r0.First();
        int st2 = r1.First();
        int st3 = r2.First();

        this->datap3 = this->data3 - st3;

        for ( int i2 = 0; i2 < bound2; ++ i2 )
        {
            int j1 = bound1 * ( i2 );  // = bound1*(i2 + j2) where j2 = 0
            this->datap3[ i2 + st3 ] = & this->data2[ j1 ] - st2;
            for ( int i1 = 0; i1 < bound1; ++ i1 )
            {
                int j0 = bound0 * ( i1 + j1 );
                this->datap3[ i2 + st3 ][ i1 + st2 ] = & this->data1[ j0 ] - st1;
            }
        }
    }

    void AllocateArray( Range r0, Range r1, Range r2, Range r3 )
    {
        int bound0 = r0.Length();
        int bound1 = r1.Length();
        int bound2 = r2.Length();
        int bound3 = r3.Length();

        int bound23   = bound2 * bound3;
        int bound123  = bound1 * bound23;
        int bound0123 = bound0 * bound123;

        this->owned4 = std::make_unique<T***[]>( static_cast<std::size_t>( bound3 ) );
        this->owned3 = std::make_unique<T**[]>( static_cast<std::size_t>( bound23 ) );
        this->owned2 = std::make_unique<T*[]>( static_cast<std::size_t>( bound123 ) );
        this->data4 = this->owned4.get();
        this->data3 = this->owned3.get();
        this->data2 = this->owned2.get();

        if ( ! shareMem )
        {
            this->owned1 = std::make_unique<T[]>( static_cast<std::size_t>( bound0123 ) );
            this->data1 = this->owned1.get();
        }
        else
        {
            this->owned1.reset();
            this->data1 = shareMem;
        }

        int st1 = r0.First();
        int st2 = r1.First();
        int st3 = r2.First();
        int st4 = r3.First();

        this->datap4 = this->data4 - st4;

        for ( int i3 = 0; i3 < bound3; ++ i3 )
        {
            int j2 = bound2 * i3;  // = bound2*(i3 + j3) where j3 = 0
            this->datap4[ i3 + st4 ] = & this->data3[ j2 ] - st3;
            for ( int i2 = 0; i2 < bound2; ++ i2 )
            {
                int j1 = bound1 * ( i2 + j2 );
                this->datap4[ i3 + st4 ][ i2 + st3 ] = & this->data2[ j1 ] - st2;
                for ( int i1 = 0; i1 < bound1; ++ i1 )
                {
                    int j0 = bound0 * ( i1 + j1 );
                    this->datap4[ i3 + st4 ][ i2 + st3 ][ i1 + st2 ] = & this->data1[ j0 ] - st1;
                }
            }
        }
    }

};

EndNameSpace

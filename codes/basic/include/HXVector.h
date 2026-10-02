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
#include <type_traits>
#include <utility>
#include <vector>

BeginNameSpace( ONEFLOW )

// HXVector is a thin adapter over std::vector used widely by the legacy code.
//
// Design notes:
// - Inherit constructors from std::vector (default, size, initializer_list, ...).
// - Do NOT use `using std::vector<T>::operator=;`. Bringing base assignment
//   into scope alongside user-declared overloads causes overload conflicts /
//   ambiguity on MSVC. Assignment is declared explicitly below instead.
// - Copy operations are constrained (C++20 requires) so move-only element
//   types such as std::unique_ptr<T> work the same way as with std::vector:
//   move is allowed, copy is not.
template < typename T >
class HXVector : public std::vector< T >
{
public:
    using std::vector< T >::vector;

    ~HXVector() = default;

    // Copy/move special members. For move-only T, copy is implicitly deleted
    // (same behavior as std::vector<T>).
    HXVector( const HXVector & ) = default;
    HXVector( HXVector && ) noexcept = default;
    HXVector & operator=( const HXVector & ) = default;
    HXVector & operator=( HXVector && ) noexcept = default;

    // Implicit conversion from std::vector: copy only when T is copyable
    // (legacy call sites depend on this conversion for copyable T).
    HXVector( const std::vector< T > & values )
        requires std::is_copy_constructible_v< T >
    : std::vector< T >( values )
    {
    }

    // Move conversion from std::vector: supports move-only T.
    HXVector( std::vector< T > && values ) noexcept
        : std::vector< T >( std::move( values ) )
    {
    }

    // Fill every element with the same value (requires copy-assignable T).
    HXVector & operator=( const T & value )
        requires std::is_copy_assignable_v< T >
    {
        for ( std::size_t i = 0; i < this->size(); ++ i )
        {
            ( *this )[ i ] = value;
        }
        return *this;
    }

    // Assign from std::vector by copy (requires copy-assignable T).
    HXVector & operator=( const std::vector< T > & values )
        requires std::is_copy_assignable_v< T >
    {
        this->assign( values.begin(), values.end() );
        return *this;
    }

    // Assign from std::vector by move (supports move-only T).
    HXVector & operator=( std::vector< T > && values ) noexcept
    {
        std::vector< T >::operator=( std::move( values ) );
        return *this;
    }
};

template < typename T >
void AllocateVector( HXVector< HXVector< T > > & data, int ni, int nj )
{
    if ( nj <= 0 ) return;
    data.resize( ni );
    for ( int i = 0; i < ni; ++ i )
    {
        data[ i ].resize( nj );
    }
}

template < typename T >
void AllocateVector( HXVector< HXVector< HXVector< T > > > & data, int ni, int nj, int nk )
{
    data.resize( ni );
    for ( int i = 0; i < ni; ++ i )
    {
        data[ i ].resize( nj );
        for ( int j = 0; j < nj; ++ j )
        {
            data[ i ][ j ].resize( nk );
        }
    }
}

template < typename T >
void Resize2D( HXVector< HXVector< T > > & data, int ni, int nj )
{
    AllocateVector( data, ni, nj );
}

EndNameSpace

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

#include <cstddef>
#include <memory>
#include <vector>

#include "Grid.h"

BeginNameSpace( ONEFLOW )

// Owning collection of grids (one entry per zone / partition piece).
using Grids = std::vector< std::unique_ptr< Grid > >;

[[nodiscard]] inline Grid & GridAt( Grids & grids, std::size_t i )
{
    return *grids[ i ];
}

[[nodiscard]] inline const Grid & GridAt( const Grids & grids, std::size_t i )
{
    return *grids[ i ];
}

[[nodiscard]] inline int GridsSize( const Grids & grids ) noexcept
{
    return static_cast< int >( grids.size() );
}

// Non-owning view of Grid pointers (caller keeps ownership).
using GridViews = std::vector< Grid * >;

[[nodiscard]] inline Grid * GridAt( GridViews & grids, std::size_t i )
{
    return grids[ i ];
}

[[nodiscard]] inline Grid * GridAt( const GridViews & grids, std::size_t i )
{
    return grids[ i ];
}

[[nodiscard]] inline int GridsSize( const GridViews & grids ) noexcept
{
    return static_cast< int >( grids.size() );
}

// Build a non-owning view from an owning collection (for read-only algorithms).
[[nodiscard]] inline GridViews AsGridViews( Grids & grids )
{
    GridViews views;
    views.reserve( grids.size() );
    for ( auto & g : grids )
    {
        views.push_back( g.get() );
    }
    return views;
}

EndNameSpace

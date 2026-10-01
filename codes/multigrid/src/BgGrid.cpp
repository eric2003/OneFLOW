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

#include "BgGrid.h"
#include "Zone.h"
#include "Mesh.h"
#include "Grid.h"
#include "StrGrid.h"
#include "UnsGrid.h"
#include "GridState.h"
#include "Multigrid.h"
#include <iostream>
#include <memory>
#include <utility>

BeginNameSpace( ONEFLOW )

namespace
{
std::unique_ptr< Grid > CloneRegistered( const char * typeName )
{
    return std::unique_ptr< Grid >( Grid::SafeClone( typeName ) );
}
}

std::unique_ptr< Grid > CreateGridUnique( int gridType )
{
    if ( gridType == ONEFLOW::UMESH )
    {
        auto grid = CloneRegistered( "UnsGrid" );
        grid->Init();
        return grid;
    }
    if ( gridType == ONEFLOW::SMESH )
    {
        auto grid = CloneRegistered( "StrGrid" );
        grid->Init();
        return grid;
    }
    std::cout << "No grid of this type\n";
    return nullptr;
}

std::unique_ptr< Grid > CreateUnsGridUnique()
{
    auto grid = std::make_unique< UnsGrid >();
    grid->Init();
    return grid;
}

std::unique_ptr< Grid > CreateStrGridUnique()
{
    auto grid = std::make_unique< StrGrid >();
    grid->Init();
    return grid;
}

Grid * CreateGrid( int gridType )
{
    return CreateGridUnique( gridType ).release();
}

Grid * CreateUnsGrid()
{
    return CreateUnsGridUnique().release();
}

Grid * CreateStrGrid()
{
    return CreateStrGridUnique().release();
}

EndNameSpace

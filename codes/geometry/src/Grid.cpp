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

#include "Grid.h"
#include "Dimension.h"
#include "NodeMesh.h"
#include "InterFace.h"
#include "SlipFace.h"
#include "DataBase.h"
#include <iostream>
#include <memory>
#include <utility>


BeginNameSpace( ONEFLOW )

namespace
{
using GridRegistry = std::map< std::string, std::unique_ptr< Grid > >;

GridRegistry & GetGridRegistry()
{
    static GridRegistry registry;
    return registry;
}
}

Grid::Grid()
{
    name = "grid";
    this->dimension = THREE_D;
    this->volBcType = -1;
}

Grid::~Grid()
{
    this->Free();
}

std::unique_ptr< Grid > Grid::SafeCloneUnique( const std::string & type )
{
    GridRegistry & registry = GetGridRegistry();
    GridRegistry::iterator iter = registry.find( type );
    if ( iter == registry.end() )
    {
        std::cout << type << " class not found" << std::endl;
        exit( 0 );
    }

    return iter->second->Clone();
}


Grid * Grid::Register( const std::string & type, std::unique_ptr< Grid > clone )
{
    GridRegistry & registry = GetGridRegistry();
    GridRegistry::iterator iter = registry.find( type );
    if ( iter != registry.end() ) return iter->second.get();

    Grid * registeredGrid = clone.get();
    registry.emplace( type, std::move( clone ) );
    return registeredGrid;
}

Grid * Grid::Register( const std::string & type, Grid * clone )
{
    return Grid::Register( type, std::unique_ptr< Grid >( clone ) );
}

void Grid::BasicInit()
{
    this->Free();
    nodeMesh  = std::make_unique< NodeMesh >();
    interFace = std::make_unique< InterFace >();
    slipFace  = std::make_unique< SlipFace >();
    dataBase  = std::make_unique< DataBase >();
}

void Grid::Free()
{
    nodeMesh.reset();
    interFace.reset();
    slipFace.reset();
    dataBase.reset();
}

void Grid::Init()
{
    this->BasicInit();
}

bool Grid::IsOneD()
{
    return this->dimension == ONEFLOW::ONE_D;
}

bool Grid::IsTwoD()
{
    return this->dimension == ONEFLOW::TWO_D;
}

bool Grid::IsThreeD()
{
    return this->dimension == ONEFLOW::THREE_D;
}

EndNameSpace

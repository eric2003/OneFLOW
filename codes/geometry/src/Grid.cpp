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

#include "Grid.h"
#include "NodeMesh.h"
#include "InterFace.h"
#include "SlipFace.h"
#include "DataBase.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

map< string, Grid * > * Grid::classMap = 0;

Grid::Grid()
{
    name = "grid";
    volBcType = -1;
    this->nodeMesh = 0;
    this->interFace = 0;
    this->dataBase = 0;
}

Grid::~Grid()
{
    this->Free();
}

Grid * Grid::SafeClone( const string & type )
{
    map < string, Grid * >::iterator iter = Grid::classMap->find( type );
    if ( iter == Grid::classMap->end() )
    {
        cout << type << " class not found" << endl;
        exit( 0 );
    }

    return iter->second->Clone();
}

Grid * Grid::Register( const string & type, Grid * clone )
{
    if ( ! Grid::classMap )
    {
        Grid::classMap = new map < string, Grid * >();
    }

    map < string, Grid * >::iterator iter = Grid::classMap->find( type );
    if ( iter == Grid::classMap->end() )
    {
        ( * Grid::classMap )[ type ] = clone;
        return clone;
    }
    else
    {
        delete clone;
        return iter->second;
    }
}

void Grid::BasicInit()
{
    nodeMesh  = new NodeMesh();
    interFace = new InterFace();
    slipFace  = new SlipFace();
    dataBase  = new DataBase();
}

void Grid::Free()
{
    delete nodeMesh;
    delete interFace;
    delete slipFace;
    delete dataBase;
}

void Grid::Init()
{
    this->BasicInit();
}

EndNameSpace
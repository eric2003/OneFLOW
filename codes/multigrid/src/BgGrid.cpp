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

#include "BgGrid.h"
#include "Zone.h"
#include "Mesh.h"
#include "Grid.h"
#include "StrGrid.h"
#include "UnsGrid.h"
#include "GridState.h"
#include "Multigrid.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

Grid * CreateGrid( int gridType )
{
    if ( gridType == ONEFLOW::UMESH )
    {
        Grid * grid = Grid::SafeClone( "UnsGrid" );
        grid->Init();
        return grid;
    }
    else if ( gridType == ONEFLOW::SMESH )
    {
        Grid * grid = Grid::SafeClone( "StrGrid" );
        grid->Init();
        return grid;
    }
    cout << "No grid of this type\n";
    return 0;
}

Grid * CreateUnsGrid()
{
    Grid * grid = new UnsGrid();
    grid->Init();
    return grid;
}

Grid * CreateStrGrid()
{
    Grid * grid = new StrGrid();
    grid->Init();
    return grid;
}

EndNameSpace
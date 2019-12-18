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

#include "GridMediator.h"
#include "Plot3D.h"
#include "Su2Grid.h"

using namespace std;

BeginNameSpace( ONEFLOW )


GridMediator::GridMediator()
{
    ;
}

GridMediator::~GridMediator()
{
    ;
}

void GridMediator::ReadGrid()
{
    if ( this->gridType == "gridgen" )
    {
        this->ReadGridgen();
    }
    else if ( this->gridType == "plot3d" )
    {
        this->ReadPlot3D();
    }
}

void GridMediator::ReadPlot3D()
{
    Plot3D::ReadPlot3D( this );
}

void GridMediator::ReadPlot3DCoor()
{
    Plot3D::ReadCoor( this );
}

void GridMediator::ReadGridgen()
{
}

EndNameSpace
/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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

#include "ScalarZone.h"
#include "MetisGrid.h"
#include "ZoneState.h"

BeginNameSpace( ONEFLOW )

int ScalarZone::nLocalZones = 0;
HXVector< ScalarGrid * > ScalarZone::scalar_grids;

ScalarZone::ScalarZone()
{
}

ScalarZone::~ScalarZone()
{
}

void ScalarZone::Allocate()
{
}

void ScalarZone::DeAllocate()
{
}

void ScalarZone::AddGrid( int zid, ScalarGrid * grid )
{
    if ( ScalarZone::scalar_grids.size() == 0 )
    {
        ScalarZone::scalar_grids.resize( ZoneState::nZones, 0 );
    }
    scalar_grids[ zid ] = grid;
}

ScalarGrid * ScalarZone::GetGrid( int iZone )
{
    return ScalarZone::scalar_grids[ iZone ];
}

ScalarGrid * ScalarZone::GetGrid()
{
    return ScalarZone::scalar_grids[ ZoneState::zid ];
}

EndNameSpace

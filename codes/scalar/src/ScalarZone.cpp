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

#include "ScalarZone.h"
#include "ScalarGrid.h"
#include "ZoneState.h"
#include <cstddef>
#include <memory>
#include <utility>


BeginNameSpace( ONEFLOW )

int ScalarZone::nLocalZones = 0;
std::vector< std::unique_ptr< ScalarGrid > > ScalarZone::scalar_grids;

ScalarZone::ScalarZone()
{
}

ScalarZone::~ScalarZone()
{
}

void ScalarZone::Allocate()
{
}

void ScalarZone::Reset()
{
    ScalarZone::scalar_grids.clear();
    ScalarZone::nLocalZones = 0;
}

void ScalarZone::AddGrid( int zid, std::unique_ptr< ScalarGrid > grid )
{
    if ( ScalarZone::scalar_grids.empty() )
    {
        ScalarZone::scalar_grids.resize( static_cast< std::size_t >( ZoneState::nZones ) );
    }
    ScalarZone::scalar_grids[ static_cast< std::size_t >( zid ) ] = std::move( grid );
}


ScalarGrid * ScalarZone::GetGrid( int iZone )
{
    return ScalarZone::scalar_grids[ static_cast< std::size_t >( iZone ) ].get();
}

ScalarGrid * ScalarZone::GetGrid()
{
    return ScalarZone::GetGrid( ZoneState::zid );
}

EndNameSpace

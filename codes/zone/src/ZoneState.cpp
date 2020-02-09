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

#include "ZoneState.h"
#include "Parallel.h"
#include "SolverDef.h"

BeginNameSpace( ONEFLOW )

int ZoneState::nLocal = 0;
int ZoneState::nZones = 0;
int ZoneState::zid = 0;
int ZoneState::szid = 0;
int ZoneState::rzid = 0;

IntField ZoneState::pid;
IntField ZoneState::zoneType;
IntField ZoneState::localZid;

ZoneState::ZoneState()
{
    ;
}

ZoneState::~ZoneState()
{
    ;
}

bool ZoneState::IsValidZone( int zoneId )
{
    return ZoneState::pid[ zoneId ] == Parallel::GetPid();
}

int ZoneState::GetZid( int iSr )
{
    if ( iSr == GREAT_SEND )
    {
        return ZoneState::szid;
    }
    return ZoneState::rzid;
}

EndNameSpace
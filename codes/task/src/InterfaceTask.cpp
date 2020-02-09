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

#include "InterfaceTask.h"
#include "Parallel.h"
#include "Zone.h"
#include "ZoneState.h"
#include "PIO.h"
#include "ActionState.h"
#include "DataBook.h"
#include "InterFace.h"

BeginNameSpace( ONEFLOW )


CUpdateInterface::CUpdateInterface()
{
    ;
}

CUpdateInterface::~CUpdateInterface()
{
    ;
}

void CUpdateInterface::Run()
{
    ActionState::dataBook = this->dataBook;
    for ( int zId = 0; zId < ZoneState::nZones; ++ zId )
    {
        ZoneState::zid = zId;

        int nNei = interFaceTopo.data[ zId ].size();

        for ( int iNei = 0; iNei < nNei; ++ iNei )
        {
            int jZone = interFaceTopo.data[ zId ][ iNei ];

            this->SwapInterfaceData( zId, jZone );
        }
    }
}

void CUpdateInterface::SwapInterfaceData( int iZone, int jZone )
{
    int sPid = ZoneState::pid[ iZone ];
    int rPid = ZoneState::pid[ jZone ];

    if ( Parallel::pid == sPid )
    {
        ZoneState::zid  = iZone;
        ZoneState::rzid = jZone;

        this->sendAction();
    }

    HXSwapData( ActionState::dataBook, sPid, rPid );

    if ( Parallel::pid == rPid )
    {
        ZoneState::zid  = jZone;
        ZoneState::szid = iZone;

        this->recvAction();
    }
}

EndNameSpace
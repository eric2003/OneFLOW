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

#include "ReadTask.h"
#include "Parallel.h"
#include "Zone.h"
#include "ZoneState.h"
#include "PIO.h"
#include "ActionState.h"
#include "DataBook.h"
#include "InterFace.h"

BeginNameSpace( ONEFLOW )

CReadFile::CReadFile()
{
}

CReadFile::~CReadFile()
{
}

void CReadFile::Run()
{
    ActionState::dataBook = this->dataBook;
    if ( Parallel::mode == 0 )
    {
        this->ServerRead();
    }
}

void CReadFile::ServerRead()
{
    fstream file;
    ActionState::file = & file;

    PIO::ParallelOpenPrj();

    for ( int zId = 0; zId < ZoneState::nZones; ++ zId )
    {
        ZoneState::zid = zId;

        this->ServerRead( this->mainAction );
    }

    PIO::ParallelClose();
}

void CReadFile::ServerRead( VoidFunc mainAction )
{
    int sPid = Parallel::serverid;
    int rPid = ZoneState::pid[ ZoneState::zid ];

    if ( Parallel::pid == sPid )
    {
        mainAction();
    }

    HXSwapData( ActionState::dataBook, sPid, rPid );

    if ( Parallel::pid == rPid )
    {
        this->action();
    }
}

EndNameSpace
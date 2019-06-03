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

#include "TaskCom.h"
#include "Parallel.h"
#include "Zone.h"
#include "ActionState.h"
#include "DataBook.h"
#include "InterFace.h"
#include <iostream>
using namespace std;

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

CWriteFile::CWriteFile()
{
}

CWriteFile::~CWriteFile()
{
}

void CWriteFile::Run()
{
    ActionState::dataBook = this->dataBook;
    if ( Parallel::mode == 0 )
    {
        this->ServerWrite();
    }
}

void CWriteFile::ServerWrite()
{
	fstream file;
    ActionState::file = & file;

	PIO::ParallelOpenPrj();

    for ( int zId = 0; zId < ZoneState::nZones; ++ zId )
    {
        ZoneState::zid = zId;

        this->ServerWrite( this->mainAction );
    }

    PIO::ParallelClose();
}

void CWriteFile::ServerWrite( VoidFunc mainAction )
{
    int sPid = ZoneState::pid[ ZoneState::zid ];
    int rPid = Parallel::serverid;

    if ( Parallel::pid == sPid )
    {
        this->action();
    }

    HXSwapData( ActionState::dataBook, sPid, rPid );

    if ( Parallel::pid == rPid )
    {
        mainAction();
    }
}

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

void Client2Server( Task * task, VoidFunc mainAction )
{
    int sPid = ZoneState::pid[ ZoneState::zid ];
    int rPid = Parallel::serverid;

    if ( Parallel::pid == sPid )
    {
        task->action();
    }

    HXSwapData( ActionState::dataBook, sPid, rPid );

	if ( Parallel::pid == rPid )
	{
        mainAction();
	}
}

void ReadBinaryFile()
{
	ActionState::dataBook->ReadFile( * ActionState::file );
}

void WriteBinaryFile()
{
	ActionState::dataBook->WriteFile( * ActionState::file );
}

void WriteAsciiFile()
{
    string str;
    ActionState::dataBook->ToString( str );
    * ActionState::file << str;
}

void WriteScreen()
{
    string str;
    ActionState::dataBook->ToString( str );
    cout << str;
}

EndNameSpace
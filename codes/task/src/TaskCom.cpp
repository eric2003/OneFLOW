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

#include "TaskCom.h"
#include "Parallel.h"
#include "Zone.h"
#include "ZoneState.h"
#include "PIO.h"
#include "ActionState.h"
#include "DataBook.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )


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
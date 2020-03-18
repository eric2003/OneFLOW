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
#include "FileMap.h"
#include "FileInfo.h"
#include "ActionState.h"
#include "DataBase.h"
#include "Task.h"
#include "TaskState.h"
#include "TaskRegister.h"
BeginNameSpace( ONEFLOW )

REGISTER_TASK( RegisterFileTask )

void RegisterFileTask()
{
    REGISTER_DATA_CLASS( SetFile );
}

void SetFile( StringField & data )
{
    string & fileNameVar = data[ 0 ];

    string fileName = GetDataValue< string >( fileNameVar );

    ios_base::openmode openMode = GetOpenMode( data[ 1 ] );

    for ( int i = 2; i < data.size(); ++ i )
    {
        openMode |= GetOpenMode( data[ i ] );
    }

    TaskState::task->fileInfo->fileName = fileName;
    TaskState::task->fileInfo->openMode = openMode;
}

ios_base::openmode GetOpenMode( const string & openModeName )
{
    if ( openModeName == "in" )
    {
        return ios_base::in;
    }
    else if ( openModeName == "out" )
    {
        return ios_base::out;
    }
    else if ( openModeName == "binary" )
    {
        return ios_base::binary;
    }
    else if ( openModeName == "trunc" )
    {
        return ios_base::trunc;
    }
    else
    {
        return ios_base::app;
    }
}

EndNameSpace
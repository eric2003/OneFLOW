/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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
#include "Iteration.h"
#include "FileUtil.h"

BeginNameSpace( ONEFLOW )

REGISTER_TASK( RegisterFileTask )

void RegisterFileTask()
{
    REGISTER_DATA_CLASS( SetFile );
}

std::string GetParallelFileName( const std::string & fileNameVar )
{
    std::string fileName = GetDataValue< std::string >( fileNameVar );

    if ( fileNameVar == "visualFile" )
    {
        int addVisualizationSteps = GetDataValue< int >( "addVisualizationSteps" );
        if ( addVisualizationSteps == 1 )
        {
            fileName = AddSymbolToFileName( fileName, Iteration::outerSteps );
        }
    }
    return fileName;
}

void SetFile( StringField & data )
{
    std::string & fileNameVar = data[ 0 ];

    std::string fileName = GetParallelFileName( fileNameVar );
    
    std::ios_base::openmode openMode = GetOpenMode( data[ 1 ] );

    for ( int i = 2; i < data.size(); ++ i )
    {
        openMode |= GetOpenMode( data[ i ] );
    }

    TaskState::task->fileInfo->fileName = fileName;
    TaskState::task->fileInfo->openMode = openMode;
}

std::ios_base::openmode GetOpenMode( const std::string & openModeName )
{
    if ( openModeName == "in" )
    {
        return std::ios_base::in;
    }
    else if ( openModeName == "out" )
    {
        return std::ios_base::out;
    }
    else if ( openModeName == "binary" )
    {
        return std::ios_base::binary;
    }
    else if ( openModeName == "trunc" )
    {
        return std::ios_base::trunc;
    }
    else
    {
        return std::ios_base::app;
    }
}

EndNameSpace

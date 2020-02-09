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

#include "LogFile.h"
#include "OStream.h"
#include "Parallel.h"
#include "Prj.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

LogFile logFile;

void OpenLogFile( int logFileIndex, fstream & file )
{
    static int ifReWrite = 0;

    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << "log/log" << logFileIndex << ".log";
    string fileName = ONEFLOW::StrIO.str();

    if ( ifReWrite == 0 )
    {
        ONEFLOW::MakePrjDir( "log" );

        ONEFLOW::OpenPrjFile( file, fileName, ios_base::out | ios_base::trunc );

        ifReWrite = 1;
    }
    else
    {
        ONEFLOW::OpenPrjFile( file, fileName, ios_base::out | ios_base::app );
    }
}

void CloseLogFile( fstream & file )
{
    file.close();
    file.clear();
}

LogFile::LogFile()
{
}

LogFile::~LogFile()
{
}

void LogFile::Open()
{
    int pid = ONEFLOW::Parallel::GetPid();
    ONEFLOW::OpenLogFile( pid, this->my_fstream );
}

void LogFile::Close()
{
    ONEFLOW::CloseLogFile( this->my_fstream );
}

EndNameSpace

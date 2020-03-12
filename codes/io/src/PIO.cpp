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

#include "PIO.h"
#include "ZoneState.h"
#include "Parallel.h"
#include "OStream.h"
#include "FileUtil.h"
#include "FileInfo.h"
#include "Stop.h"
#include "Prj.h"
#include "ActionState.h"
#include "Task.h"
#include "TaskState.h"
#include <iostream>
#include <string>
using namespace std;

BeginNameSpace( ONEFLOW )


PIO::PIO()
{
    ;
}

PIO::~PIO()
{
    ;
}

string PIO::GetDirName( const string & fileName )
{
    size_t pos = fileName.find_last_of("\\/");
    if ( string::npos == pos )
    {
        return "";
    }
    else
    {
        return fileName.substr(0, pos);
    }
}

void PIO::ParallelOpen( fstream & file, const string & fileName, const ios_base::openmode & openMode )
{
    if ( Parallel::pid != Parallel::GetFid() ) return;

    PIO::Open( file, fileName, openMode );
}

void PIO::ParallelOpenPrj( fstream & file, const string & fileName, const ios_base::openmode & openMode )
{
    if ( Parallel::pid != Parallel::GetFid() ) return;

    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << PrjStatus::prjBaseDir << fileName;

    string prjFileName = ONEFLOW::StrIO.str();
    string prj_dir = PIO::GetDirName( prjFileName );

    if ( ! DirExist( prj_dir ) )
    {
        MakeDir( prj_dir );
    }
    PIO::Open( file, prjFileName, openMode );
}

void PIO::ParallelOpenPrj()
{
    PIO::ParallelOpenPrj( * ActionState::file, TaskState::task->fileInfo->fileName, TaskState::task->fileInfo->openMode );
}

void PIO::ParallelClose()
{
    PIO::ParallelClose( * ActionState::file );
}

void PIO::ParallelClose( fstream & file )
{
    if ( Parallel::pid != Parallel::GetFid() ) return;

    PIO::Close( file );
}

void PIO::Open( fstream & file, const string & fileName, const ios_base::openmode & openMode )
{
    file.open( fileName.c_str(), openMode );
    if ( ! file )
    {
        cout << "could not open " << fileName << endl;
        Stop("");
    }
}

void PIO::Close( fstream & file )
{
    file.close();
    file.clear();
}


EndNameSpace
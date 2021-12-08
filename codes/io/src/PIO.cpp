/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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

BeginNameSpace( ONEFLOW )


PIO::PIO()
{
    ;
}

PIO::~PIO()
{
    ;
}

void PIO::ParallelOpen( std::fstream & file, const std::string & fileName, const std::ios_base::openmode & openMode )
{
    if ( Parallel::pid != Parallel::GetFid() ) return;

    Prj::OpenFile( file, fileName, openMode );
}

void PIO::OpenPrjFile( std::fstream & file, const std::string & fileName, const std::ios_base::openmode & openMode )
{
    if ( Parallel::pid != Parallel::GetFid() ) return;

    Prj::OpenPrjFile( file, fileName, openMode );
}

void PIO::OpenPrjFile()
{
    PIO::OpenPrjFile( * ActionState::file, TaskState::task->fileInfo->fileName, TaskState::task->fileInfo->openMode );
}

void PIO::CloseFile()
{
    PIO::CloseFile( * ActionState::file );
}

void PIO::CloseFile( std::fstream & file )
{
    if ( Parallel::pid != Parallel::GetFid() ) return;

    Prj::CloseFile( file );
}


EndNameSpace

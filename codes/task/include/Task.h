/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
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
#pragma once

#pragma once

#include "HXDefine.h"

#include <memory>
#include <string>

BeginNameSpace( ONEFLOW )

using TaskFunction = void ( * )( void );

class DataBook;
class FileInfo;

class Task
{
public:
    Task();

    virtual ~Task();

    Task( const Task & ) = delete;
    Task & operator=( const Task & ) = delete;

    Task( Task && ) noexcept;
    Task & operator=( Task && ) noexcept;

public:
    virtual void Run() {}

public:
    int taskId = -1;
    std::string taskName;

    TaskFunction action = nullptr;
    TaskFunction sendAction = nullptr;
    TaskFunction recvAction = nullptr;

    std::unique_ptr< DataBook > dataBook;
    std::unique_ptr< FileInfo > fileInfo;
};

EndNameSpace

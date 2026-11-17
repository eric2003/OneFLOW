/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software; you can redistribute it and/or modify it
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

#include "NamespaceMacros.h"

#include <string>
#include <vector>

BeginNameSpace( ONEFLOW )

struct ProcessResult
{
    bool started{ false };
    int exitCode{ -1 };
    std::string output;
    std::string errorMessage;
};

// Starts an executable directly without invoking a command shell.
// Arguments are UTF-8 strings and are passed as distinct process arguments.
// Standard output and standard error are captured into one diagnostic stream.
ProcessResult RunProcess(
    const std::string & executable,
    const std::vector< std::string > & arguments = {},
    const std::string & workingDirectory = {} );

EndNameSpace

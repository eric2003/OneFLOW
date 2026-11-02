/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
    This file is part of OneFLOW.

    OneFLOW is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OneFLOW is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.
----------------------------------------------------------------------------*/

#include "Prj.h"
#include "SimulationCommandLine.h"

#include <stdexcept>
#include <vector>

BeginNameSpace( ONEFLOW )

CmdLineOptions Prj::ParseCmdLineArgs(
    const std::vector<std::string> & args )
{
    const SimulationCommandLineOptions parsed =
        ParseSimulationCommandLineOptions( args );

    CmdLineOptions options;
    options.debug = parsed.debug;
    options.caseDir = parsed.caseDirectories.front();
    options.caseDirs = parsed.caseDirectories;
    return options;
}

EndNameSpace

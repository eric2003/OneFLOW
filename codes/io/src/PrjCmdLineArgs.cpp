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

// Pure command-line parsing for Prj::ProcessCmdLineArgs.
// Deliberately kept dependency-free (no Fatal/OStream/FileUtils) so it can
// be compiled and unit-tested in isolation from the rest of Prj.cpp.

#include "Prj.h"
#include <stdexcept>

BeginNameSpace( ONEFLOW )

CmdLineOptions Prj::ParseCmdLineArgs( const std::vector<std::string> & args )
{
    if ( args.size() < 3 )
    {
        throw std::invalid_argument(
            "Prj::ParseCmdLineArgs: expected at least 3 arguments "
            "(exe, mode, prjName), got " + std::to_string( args.size() ) );
    }

    CmdLineOptions opt;
    opt.debug   = ( args[ 1 ] == "d" );
    opt.prjName = args[ 2 ];
    return opt;
}

EndNameSpace
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

// Adapts OneFLOW positional arguments to the reusable CLI parsing core.

#include "Prj.h"
#include "CliParser.h"

#include <stdexcept>
#include <vector>

BeginNameSpace( ONEFLOW )

CmdLineOptions Prj::ParseCmdLineArgs(
    const std::vector<std::string> & args )
{
    CliParser parser( "OneFLOW" );
    parser.AddPositional(
        "mode", "Run mode (d selects debug)", 1, 1 );
    parser.AddPositional(
        "caseDir", "Case directory; repeated values run repeatedly",
        1, CliParser::Unlimited );

    if ( args.empty() )
    {
        throw std::invalid_argument(
            "Missing executable name\n" + parser.Usage() );
    }

    const std::vector<std::string> positionalArgs(
        args.begin() + 1, args.end() );
    const ParsedCommandLine parsed = parser.Parse( positionalArgs );

    CmdLineOptions options;
    options.debug = ( parsed.GetOne( "mode" ) == "d" );
    options.caseDir = parsed.GetAll( "caseDir" ).front();
    options.caseDirs = parsed.GetAll( "caseDir" );
    return options;
}

EndNameSpace

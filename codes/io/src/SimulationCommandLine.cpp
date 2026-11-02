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

    OneFLOW is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.
----------------------------------------------------------------------------*/

#include "CliParser.h"
#include "NamespaceMacros.h"
#include "SimulationCommandLine.h"

#include <stdexcept>
#include <vector>

BeginNameSpace( ONEFLOW )

SimulationInvocationKind ClassifySimulationInvocation(
    const std::vector<std::string> & args ) noexcept
{
    if ( args.size() == 1 )
    {
        return SimulationInvocationKind::LightweightTest;
    }
    if ( args.size() == 2 )
    {
        return SimulationInvocationKind::InvalidArgumentCount;
    }
    return SimulationInvocationKind::FullSimulation;
}

SimulationCommandLineOptions ParseSimulationCommandLineOptions(
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

    SimulationCommandLineOptions options;
    options.debug = ( parsed.GetOne( "mode" ) == "d" );
    options.caseDirectories = parsed.GetAll( "caseDir" );
    return options;
}

EndNameSpace
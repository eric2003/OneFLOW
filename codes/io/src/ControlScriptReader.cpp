/*---------------------------------------------------------------------------*\\
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

\\*---------------------------------------------------------------------------*/

#include "ParaFile.h"
#include "Prj.h"

#include <iostream>
#include <string>

namespace ONEFLOW {
namespace {

void ReadDiscoveredCaseScripts( const std::string & caseDir )
{
    const ScriptFileDiscovery discovery = DiscoverCaseScriptFiles( caseDir );
    for ( const std::string & diagnostic : discovery.diagnostics )
        std::cerr << "[OneFLOW parameter files] " << diagnostic << std::endl;
    ReadMultiScriptFiles( discovery.filePaths );
}

} // namespace

void ReadPrjScript()
{
    ReadDiscoveredCaseScripts( Prj::prjBaseDir );
}

void ReadPrjScript( const std::string & caseDir )
{
    ReadDiscoveredCaseScripts( caseDir );
}

} // namespace ONEFLOW

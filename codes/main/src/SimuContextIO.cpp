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
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.
----------------------------------------------------------------------------*/

#include "SimuContextIO.h"

#include "Prj.h"

#include <filesystem>
#include <string>

BeginNameSpace( ONEFLOW )

namespace SimuContextDetail
{

void ApplyCaseIOSettings( const CaseIOSettings & settings )
{
    if ( ! settings.projectRoot.empty() )
    {
        Prj::SetProjectBaseDir( settings.projectRoot );
    }

    Prj::SetIOBaseDirs( settings.inputRoot, settings.outputRoot );
    Prj::SetDomainIOBaseDirs(
        settings.meshInputRoot,
        settings.flowInputRoot,
        settings.meshOutputRoot,
        settings.flowOutputRoot );

    const bool outputsUseCaseRoot =
        Prj::AreSameDirectoryPath(
            Prj::GetMeshOutputBaseDir(), Prj::prjBaseDir )
        && Prj::AreSameDirectoryPath(
            Prj::GetFlowOutputBaseDir(), Prj::prjBaseDir );
    std::string ioPathInfoFile = outputsUseCaseRoot
        ? "results/io_paths.txt"
        : "io_paths.txt";
    if ( ! settings.ioPathInfoFileName.empty() )
    {
        ioPathInfoFile = settings.ioPathInfoFileName;
    }

    const std::string meshInfoPath =
        Prj::GetMeshOutputFileName( ioPathInfoFile );
    Prj::WriteIOPathInfo( meshInfoPath );

    const std::string flowInfoPath =
        Prj::GetFlowOutputFileName( ioPathInfoFile );
    if ( std::filesystem::path( flowInfoPath ).lexically_normal()
        != std::filesystem::path( meshInfoPath ).lexically_normal() )
    {
        Prj::WriteIOPathInfo( flowInfoPath );
    }
}

}

EndNameSpace

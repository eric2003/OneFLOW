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

#pragma once

#include "DataBaseType.h"
#include "DataObject.h"
#include "DataPara.h"

#include <stdexcept>
#include <string>

BeginNameSpace( ONEFLOW )

namespace SimuContextDetail
{

// Snapshot of case I/O settings. Read and validate the complete snapshot before
// applying any of its values to process-global Prj state.
struct CaseIOSettings
{
    std::string projectRoot;
    std::string inputRoot;
    std::string outputRoot;
    std::string meshInputRoot;
    std::string flowInputRoot;
    std::string meshOutputRoot;
    std::string flowOutputRoot;
    std::string ioPathInfoFileName;
};

inline bool ShouldUseLegacySourceGridCaseDir( const std::string & meshInputRoot )
{
    return meshInputRoot.empty() || meshInputRoot == "default";
}

inline std::string ApplyLegacySourceGridCaseDir(
    const std::string & meshInputRoot,
    const std::string & sourceGridCaseDir )
{
    if ( ! ShouldUseLegacySourceGridCaseDir( meshInputRoot ) )
    {
        return meshInputRoot;
    }

    if ( sourceGridCaseDir.empty() || sourceGridCaseDir == "default" )
    {
        return meshInputRoot;
    }

    return sourceGridCaseDir;
}

inline std::string ReadOptionalStringParameter(
    const DataPara & parameters, const std::string & name )
{
    const DataEntry * entry = parameters.FindDataEntry( name );

    if ( entry == nullptr )
    {
        return "";
    }

    if ( entry->GetType() != HX_STRING )
    {
        throw std::invalid_argument(
            "Control parameter \"" + name + "\" must be a string" );
    }

    return GetDataValue< std::string >( &entry->GetDataObject(), 0 );
}

// Missing optional keys resolve to empty strings. Present keys must have
// string type. Read all setup-time settings here so invalid configuration is
// rejected before SetupCaseEnvironment applies configured roots.
inline CaseIOSettings ReadCaseIOSettings( const DataPara & parameters )
{
    CaseIOSettings settings;
    settings.projectRoot = ReadOptionalStringParameter( parameters, "projectRoot" );
    settings.inputRoot = ReadOptionalStringParameter( parameters, "inputRoot" );
    settings.outputRoot = ReadOptionalStringParameter( parameters, "outputRoot" );
    settings.meshInputRoot = ReadOptionalStringParameter( parameters, "meshInputRoot" );
    settings.flowInputRoot = ReadOptionalStringParameter( parameters, "flowInputRoot" );
    settings.meshOutputRoot = ReadOptionalStringParameter( parameters, "meshOutputRoot" );
    settings.flowOutputRoot = ReadOptionalStringParameter( parameters, "flowOutputRoot" );
    settings.ioPathInfoFileName =
        ReadOptionalStringParameter( parameters, "ioPathInfoFileName" );

    if ( ShouldUseLegacySourceGridCaseDir( settings.meshInputRoot ) )
    {
        const std::string sourceGridCaseDir =
            ReadOptionalStringParameter( parameters, "sourceGridCaseDir" );
        settings.meshInputRoot = ApplyLegacySourceGridCaseDir(
            settings.meshInputRoot, sourceGridCaseDir );
    }

    return settings;
}

}

EndNameSpace

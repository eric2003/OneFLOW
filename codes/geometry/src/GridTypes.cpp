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

#include "GridTypes.h"
#include "DataBase.h"
#include <exception>

BeginNameSpace( ONEFLOW )

std::string GridConfig::ResolveSourceCaseDir( const std::string & currentCaseDir ) const
{
    return sourceCaseDir.empty() ? currentCaseDir : sourceCaseDir;
}

std::string GridConfig::GetSourceCaseDir()
{
    try
    {
        return GetDataValue< std::string >( "sourceGridCaseDir" );
    }
    catch ( const std::exception & )
    {
        return std::string();
    }
}

GridConfig GridConfig::FromDataBase()
{
    GridConfig cfg;

    cfg.sourceFile = GetDataValue< std::string >( "sourceGridFileName" );

    // Keep the default empty so existing cases continue to read their own grid.
    // The explicit source case is consumed later by the runtime grid reader.
    cfg.sourceCaseDir = GridConfig::GetSourceCaseDir();
    cfg.layoutFile = GetDataValue< std::string >( "gridLayoutFileName" );

    cfg.bcFile         = GetDataValue< std::string >( "sourceGridBcName" );
    cfg.targetFile     = GetDataValue< std::string >( "targetGridFileName" );

    cfg.sourceType = ParseGridFileType( GetDataValue< std::string >( "sourceGridType" ) );
    cfg.targetType = ParseGridFileType( GetDataValue< std::string >( "targetGridType" ) );
    cfg.topology = ParseGridTopology( GetDataValue< std::string >( "topoType" ) );

    cfg.multiBlock = GetDataValue< int >( "multiBlock" ) != 0;
    cfg.axisDirection = GetDataValue< int >( "axis_dir" ) == 1
        ? GridAxisDirection::ZToY
        : GridAxisDirection::Y;

    // Partition-only parameters are required only when the partition workflow is selected.
    // Other grid workflows may legitimately omit these database entries.
    const int rawObj = GetDataValue< int >( "gridObj" );
    if ( auto parsed = ParseGridObjective( rawObj ) )
    {
        cfg.objective = *parsed;
    }
    else
    {
        cfg.objective = static_cast< GridObjective >( rawObj );
    }

    if ( cfg.objective == GridObjective::Partition )
    {
        cfg.partitionFile = GetDataValue< std::string >( "part_uns_file" );
        cfg.partitionType = GetDataValue< int >( "partition_type" );
    }

    try
    {
        cfg.ignoreNoBoundary = GetDataValue< int >( "ignoreNoBc" ) != 0;
    }
    catch ( const std::exception & )
    {
        // Keep the default when legacy boundary control is not configured.
    }

    cfg.scale = GetDataValue< Real >( "gridScale" );

    const int generationId = GetDataValue< int >( "igene" );
    cfg.generationType = ParseGridGenerationType( generationId );

    // Prefer GetDataPointer over CopyArray so this TU only needs DataBase.h.
    cfg.translate = { 0.0, 0.0, 0.0 };
    try
    {
        Real * p = GetDataPointer< Real >( "gridTrans" );
        if ( p )
        {
            for ( size_t i = 0; i < 3; ++i )
            {
                cfg.translate[ i ] = p[ i ];
            }
        }
    }
    catch ( const std::exception & )
    {
        // Missing gridTrans: keep zeros.
    }

    return cfg;
}

EndNameSpace

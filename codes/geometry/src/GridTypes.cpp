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

    cfg.bcFile     = GetDataValue< std::string >( "sourceGridBcName" );
    cfg.targetFile = GetDataValue< std::string >( "targetGridFileName" );

    // These values belong to specific grid workflows. Keep their defaults when
    // a smaller workflow does not register the corresponding database entries.
    try
    {
        cfg.sourceType =
            ParseGridFileType( GetDataValue< std::string >( "sourceGridType" ) );
    }
    catch ( const std::exception & )
    {
    }

    try
    {
        cfg.targetType =
            ParseGridFileType( GetDataValue< std::string >( "targetGridType" ) );
    }
    catch ( const std::exception & )
    {
    }

    try
    {
        cfg.topology =
            ParseGridTopology( GetDataValue< std::string >( "topoType" ) );
    }
    catch ( const std::exception & )
    {
    }

    try
    {
        cfg.assemblyMode = GetDataValue< int >( "multiBlock" ) != 0
            ? GridAssemblyMode::PerZone
            : GridAssemblyMode::AggregateZones;
    }
    catch ( const std::exception & )
    {
    }

    try
    {
        cfg.axisDirection = GetDataValue< int >( "axis_dir" ) == 1
            ? GridAxisDirection::ZToY
            : GridAxisDirection::Y;
    }
    catch ( const std::exception & )
    {
    }

    // Objective is optional because generation and output workflows do not need
    // to select a conversion or partition objective.
    try
    {
        const int rawObj = GetDataValue< int >( "gridObj" );
        if ( auto parsed = ParseGridObjective( rawObj ) )
        {
            cfg.objective = *parsed;
        }
        else
        {
            cfg.objective = static_cast< GridObjective >( rawObj );
        }
    }
    catch ( const std::exception & )
    {
    }

    // Partition-only parameters are required only when the partition workflow is selected.
    // Other grid workflows may legitimately omit these database entries.
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

    try
    {
        cfg.scale = GetDataValue< Real >( "gridScale" );
    }
    catch ( const std::exception & )
    {
        // Keep the unit scale for workflows that do not transform a source grid.
    }

    try
    {
        cfg.generationId = GetDataValue< int >( "igene" );
    }
    catch ( const std::exception & )
    {
        // Keep generationType empty when classic-grid generation is not selected.
    }

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

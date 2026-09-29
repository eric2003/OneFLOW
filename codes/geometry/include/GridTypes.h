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

#pragma once

#include "HXDefine.h"
#include <array>
#include <cctype>
#include <optional>
#include <string>
#include <string_view>

BeginNameSpace( ONEFLOW )

// ---------------------------------------------------------------------------
// Offline grid factory objectives (replaces magic int gridObj).
// Values match the historical DataBase "gridObj" integers for compatibility.
// ---------------------------------------------------------------------------
enum class GridObjective : int
{
    GenerateClassic = 0,  // classic shapes then convert
    ConvertOnly     = 1,  // format conversion only
    GenerateInp     = 2,  // multi-block structured input
    Partition       = 3   // domain decomposition
};

// ---------------------------------------------------------------------------
// Grid file formats used by ConvertGrid pipelines.
// ---------------------------------------------------------------------------
enum class GridFileType
{
    Plot3D,
    SU2,
    CGNS,
    OneFLOW,
    Gridgen,
    Unknown
};

// ---------------------------------------------------------------------------
// Runtime grid operations registered under system/grid (taskMap / funcMap).
// Token strings match historical map files (case-insensitive parse).
// ---------------------------------------------------------------------------
enum class GridOp
{
    CalcMetrics,
    SwapCellCenter,
    FillWallStruct,
    CalcWallDist,
    ReadWallDist,
    WriteWallDist,
    AllocWallDist
};

// Canonical token table (order matches enum underlying values when sequential).
// Use for validation / documentation of system/grid maps.
inline constexpr std::array< std::string_view, 7 > kGridOpTokens = {
    "CALC_METRICS",
    "SWAP_CELLCENTER",
    "FILL_WALL_STRUCT",
    "CALC_WALL_DIST",
    "READ_WALL_DIST",
    "WRITE_WALL_DIST",
    "ALLOCATE_WALL_DIST"
};

// ---------------------------------------------------------------------------
// Immutable-ish configuration snapshot loaded once from DataBase.
// ---------------------------------------------------------------------------
struct GridConfig
{
    GridObjective objective{ GridObjective::ConvertOnly };
    GridFileType  sourceType{ GridFileType::Unknown };
    GridFileType  targetType{ GridFileType::Unknown };
    std::string   sourceFile;
    // Empty means the current case; otherwise this identifies the case that owns the source grid.
    std::string   sourceCaseDir;
    std::string   bcFile;
    std::string   targetFile;
    std::string   topo;
    int           multiBlock{ 0 };
    int           axisDir{ 0 };
    Real          scale{ 1.0 };
    std::array< Real, 3 > translate{};

    // Implemented in GridTypes.cpp (needs DataBase).
    static GridConfig FromDataBase();

    // Return the optional external source case for the current grid configuration.
    static std::string GetSourceCaseDir();
};

// ---------------------------------------------------------------------------
// Parsing / formatting helpers - header-only to avoid extra link deps.
// ---------------------------------------------------------------------------

[[nodiscard]] constexpr std::optional< GridObjective >
ParseGridObjective( int value ) noexcept
{
    switch ( value )
    {
        case 0: return GridObjective::GenerateClassic;
        case 1: return GridObjective::ConvertOnly;
        case 2: return GridObjective::GenerateInp;
        case 3: return GridObjective::Partition;
        default: return std::nullopt;
    }
}

namespace grid_types_detail
{
inline bool EqualIgnoreCase( std::string_view a, std::string_view b ) noexcept
{
    if ( a.size() != b.size() ) return false;
    for ( size_t i = 0; i < a.size(); ++i )
    {
        if ( std::tolower( static_cast< unsigned char >( a[ i ] ) ) !=
             std::tolower( static_cast< unsigned char >( b[ i ] ) ) )
        {
            return false;
        }
    }
    return true;
}
} // namespace grid_types_detail

[[nodiscard]] inline GridFileType ParseGridFileType( std::string_view name ) noexcept
{
    if ( grid_types_detail::EqualIgnoreCase( name, "plot3d" ) )  return GridFileType::Plot3D;
    if ( grid_types_detail::EqualIgnoreCase( name, "su2" ) )     return GridFileType::SU2;
    if ( grid_types_detail::EqualIgnoreCase( name, "cgns" ) )    return GridFileType::CGNS;
    if ( grid_types_detail::EqualIgnoreCase( name, "oneflow" ) ) return GridFileType::OneFLOW;
    if ( grid_types_detail::EqualIgnoreCase( name, "gridgen" ) ) return GridFileType::Gridgen;
    return GridFileType::Unknown;
}

[[nodiscard]] inline std::string_view ToString( GridFileType type ) noexcept
{
    switch ( type )
    {
        case GridFileType::Plot3D:  return "plot3d";
        case GridFileType::SU2:     return "su2";
        case GridFileType::CGNS:    return "cgns";
        case GridFileType::OneFLOW: return "oneflow";
        case GridFileType::Gridgen: return "gridgen";
        default:                    return "unknown";
    }
}

[[nodiscard]] inline std::string_view ToString( GridObjective obj ) noexcept
{
    switch ( obj )
    {
        case GridObjective::GenerateClassic: return "GenerateClassic";
        case GridObjective::ConvertOnly:     return "ConvertOnly";
        case GridObjective::GenerateInp:     return "GenerateInp";
        case GridObjective::Partition:       return "Partition";
        default:                             return "Unknown";
    }
}

[[nodiscard]] inline std::string_view ToString( GridOp op ) noexcept
{
    switch ( op )
    {
        case GridOp::CalcMetrics:    return kGridOpTokens[ 0 ];
        case GridOp::SwapCellCenter: return kGridOpTokens[ 1 ];
        case GridOp::FillWallStruct: return kGridOpTokens[ 2 ];
        case GridOp::CalcWallDist:   return kGridOpTokens[ 3 ];
        case GridOp::ReadWallDist:   return kGridOpTokens[ 4 ];
        case GridOp::WriteWallDist:  return kGridOpTokens[ 5 ];
        case GridOp::AllocWallDist:  return kGridOpTokens[ 6 ];
        default:                     return "UNKNOWN_GRID_OP";
    }
}

[[nodiscard]] inline std::optional< GridOp > ParseGridOp( std::string_view name ) noexcept
{
    if ( grid_types_detail::EqualIgnoreCase( name, kGridOpTokens[ 0 ] ) ) return GridOp::CalcMetrics;
    if ( grid_types_detail::EqualIgnoreCase( name, kGridOpTokens[ 1 ] ) ) return GridOp::SwapCellCenter;
    if ( grid_types_detail::EqualIgnoreCase( name, kGridOpTokens[ 2 ] ) ) return GridOp::FillWallStruct;
    if ( grid_types_detail::EqualIgnoreCase( name, kGridOpTokens[ 3 ] ) ) return GridOp::CalcWallDist;
    if ( grid_types_detail::EqualIgnoreCase( name, kGridOpTokens[ 4 ] ) ) return GridOp::ReadWallDist;
    if ( grid_types_detail::EqualIgnoreCase( name, kGridOpTokens[ 5 ] ) ) return GridOp::WriteWallDist;
    if ( grid_types_detail::EqualIgnoreCase( name, kGridOpTokens[ 6 ] ) ) return GridOp::AllocWallDist;
    return std::nullopt;
}

// Returns true if token is a known grid operation (case-insensitive).
[[nodiscard]] inline bool IsKnownGridOpToken( std::string_view name ) noexcept
{
    return ParseGridOp( name ).has_value();
}

EndNameSpace

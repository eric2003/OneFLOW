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

#include "GridConversion.h"
#include "CgnsFactory.h"
#include "GridMediator.h"
#include "Su2Grid.h"
#include "Plot3D.h"
#include "Boundary.h"
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

BeginNameSpace( ONEFLOW )

namespace
{
    using Converter = void ( * )(
        const GridConfig &,
        const std::string & );

    struct ConverterEntry
    {
        GridFileType sourceType;
        Converter convert;
    };


    void ConvertGmshToCgnsPlaceholder(
        const GridConfig & config,
        const std::string & caseDir )
    {
        const std::filesystem::path sourceRoot = config.ResolveSourceCaseDir( caseDir );
        const std::filesystem::path sourceName( config.sourceFile );
        const std::filesystem::path source = sourceName.is_absolute()
            ? sourceName : sourceRoot / sourceName;
        const std::filesystem::path targetName( config.targetFile );
        const std::filesystem::path target = targetName.is_absolute()
            ? targetName : std::filesystem::path( caseDir ) / targetName;

        std::ifstream input( source );
        if ( ! input ) throw std::runtime_error( "Cannot open Gmsh input: " + source.string() );

        long long nodeCount = 0, elementCount = 0;
        int versionMajor = 0, versionMinor = 0;
        bool hasFormat = false, hasNodes = false, hasElements = false;
        std::string line;
        while ( std::getline( input, line ) )
        {
            // Gmsh files may use CRLF even when read on Linux.
            if ( ! line.empty() && line.back() == '\r' ) line.pop_back();
            if ( line == "$MeshFormat" )
            {
                if ( ! std::getline( input, line ) ) throw std::runtime_error( "Truncated Gmsh format header" );
                std::istringstream row( line );
                double version = 0.0;
                int fileType = -1, dataSize = 0;
                if ( ! ( row >> version >> fileType >> dataSize ) ||
                     fileType != 0 || version < 2.0 || version >= 3.0 )
                    throw std::runtime_error( "Initial Gmsh reader supports ASCII MSH 2.x only" );
                versionMajor = 2;
                versionMinor = static_cast< int >( ( version - 2.0 ) * 10.0 + 0.5 );
                hasFormat = true;
            }
            else if ( line == "$Nodes" )
            {
                if ( ! std::getline( input, line ) ) throw std::runtime_error( "Truncated Gmsh nodes section" );
                std::istringstream count( line );
                if ( ! ( count >> nodeCount ) || nodeCount < 1 ) throw std::runtime_error( "Invalid Gmsh node count" );
                for ( long long i = 0; i < nodeCount; ++i )
                {
                    if ( ! std::getline( input, line ) ) throw std::runtime_error( "Truncated Gmsh node records" );
                    std::istringstream node( line );
                    long long id = 0;
                    double x = 0.0, y = 0.0, z = 0.0;
                    if ( ! ( node >> id >> x >> y >> z ) ) throw std::runtime_error( "Invalid Gmsh node record" );
                }
                hasNodes = true;
            }
            else if ( line == "$Elements" )
            {
                if ( ! std::getline( input, line ) ) throw std::runtime_error( "Truncated Gmsh elements section" );
                std::istringstream count( line );
                if ( ! ( count >> elementCount ) || elementCount < 1 ) throw std::runtime_error( "Invalid Gmsh element count" );
                for ( long long i = 0; i < elementCount; ++i )
                {
                    if ( ! std::getline( input, line ) ) throw std::runtime_error( "Truncated Gmsh element records" );
                    std::istringstream element( line );
                    long long id = 0, ignored = 0;
                    int type = 0, tags = 0;
                    if ( ! ( element >> id >> type >> tags ) || tags < 0 ) throw std::runtime_error( "Invalid Gmsh element record" );
                    for ( int tag = 0; tag < tags; ++tag )
                        if ( ! ( element >> ignored ) ) throw std::runtime_error( "Truncated Gmsh element tags" );
                    int nodes = 0;
                    while ( element >> ignored ) ++nodes;
                    if ( nodes == 0 ) throw std::runtime_error( "Gmsh element has no node references" );
                }
                hasElements = true;
            }
        }
        if ( ! hasFormat || ! hasNodes || ! hasElements )
            throw std::runtime_error( "Gmsh input requires $MeshFormat, $Nodes and $Elements sections" );

        if ( ! target.parent_path().empty() ) std::filesystem::create_directories( target.parent_path() );
        std::ofstream output( target, std::ios::trunc );
        if ( ! output ) throw std::runtime_error( "Cannot create placeholder output: " + target.string() );
        // This deliberately marked text artifact is NOT a standards-compliant CGNS file.
        output << "ONEFLOW_CGNS_PLACEHOLDER 1\n"
               << "status=NOT_A_CGNS_FILE\n"
               << "source_format=gmsh_msh_2_ascii\n"
               << "gmsh_version=" << versionMajor << "." << versionMinor << "\n"
               << "node_count=" << nodeCount << "\n"
               << "element_count=" << elementCount << "\n";
        if ( ! output ) throw std::runtime_error( "Failed writing placeholder output: " + target.string() );
    }

    void ConvertGmsh( const GridConfig & config, const std::string & caseDir )
    {
        if ( config.targetType != GridFileType::CGNS )
            throw std::invalid_argument( "Initial Gmsh conversion supports targetGridType = cgns only" );
        ConvertGmshToCgnsPlaceholder( config, caseDir );
    }

    void ConvertPlot3DToOneFLOW(
        const GridConfig & config,
        const std::string & /*caseDir*/ )
    {
        CgnsFactory cgnsFactory;
        cgnsFactory.CommonToOneFlowGrid( config );
    }

    void ConvertPlot3DToCGNS(
        const GridConfig & config,
        const std::string & caseDir )
    {
        CgnsFactory cgnsFactory;
        GridMediatorGroup gridMediatorGroup;
        Plot3D::Plot3DToCgns( &gridMediatorGroup, config, caseDir );
        cgnsFactory.DumpCgnsGrid( gridMediatorGroup );
    }

    void ConvertPlot3D(
        const GridConfig & config,
        const std::string & caseDir )
    {
        switch ( config.targetType )
        {
            case GridFileType::OneFLOW:
                ConvertPlot3DToOneFLOW( config, caseDir );
                return;
            case GridFileType::CGNS:
                ConvertPlot3DToCGNS( config, caseDir );
                return;
            default:
                throw std::invalid_argument(
                    std::string( "Unsupported Plot3D target type: " ) +
                    std::string( ToString( config.targetType ) ) );
        }
    }

    void ConvertSU2(
        const GridConfig & config,
        const std::string & caseDir )
    {
        Su2Grid su2Grid;
        su2Grid.Su2ToOneFlowGrid( config, caseDir );
    }

    void ConvertCGNS(
        const GridConfig & config,
        const std::string & caseDir )
    {
        CgnsFactory cgnsFactory;
        cgnsFactory.GenerateGrid( config, caseDir );
    }

    constexpr ConverterEntry kConverters[] = {
        { GridFileType::Plot3D, &ConvertPlot3D },
        { GridFileType::SU2,    &ConvertSU2 },
        { GridFileType::CGNS,   &ConvertCGNS },
        { GridFileType::Gmsh,   &ConvertGmsh },
    };
}

namespace
{
bool IsSupportedConversion( GridFileType sourceType, GridFileType targetType )
{
    switch ( sourceType )
    {
        case GridFileType::Plot3D:
            return targetType == GridFileType::OneFLOW ||
                   targetType == GridFileType::CGNS;
        case GridFileType::SU2:
            return targetType == GridFileType::OneFLOW;
        case GridFileType::CGNS:
            return targetType == GridFileType::OneFLOW ||
                   targetType == GridFileType::CGNS;
        case GridFileType::Gmsh:
            return targetType == GridFileType::CGNS;
        default:
            return false;
    }
}
} // namespace

void ConvertGrid(
    const GridConfig & config,
    const std::string & caseDir )
{
    if ( ! IsSupportedConversion( config.sourceType, config.targetType ) )
    {
        throw std::invalid_argument(
            std::string( "Unsupported grid conversion: " ) +
            std::string( ToString( config.sourceType ) ) + " -> " +
            std::string( ToString( config.targetType ) ) );
    }
    // Region IDs are local to one conversion and must not inherit names from
    // a previous case executed in the same process.
    RegionNameMap::Reset();

    for ( const auto & entry : kConverters )
    {
        if ( entry.sourceType == config.sourceType )
        {
            entry.convert( config, caseDir );
            return;
        }
    }

    throw std::invalid_argument(
        std::string( "Unsupported source grid type: " ) +
        std::string( ToString( config.sourceType ) ) );
}

EndNameSpace

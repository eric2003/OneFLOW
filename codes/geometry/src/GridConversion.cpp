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
#include "GmshMeshReader.h"
#include "StlMeshIO.h"
#include "ProcessRunner.h"
#include "PointLocator.h"
#include "CgnsFactory.h"
#include "CgnsZbase.h"
#include "CgnsBase.h"
#include "CgnsZone.h"
#include "GridMediator.h"
#include "Su2Grid.h"
#include "UnstructuredMeshData.h"
#include "UnstructuredCgnsZoneWriter.h"
#include "Plot3D.h"
#include "Boundary.h"
#include <filesystem>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iterator>
#include <sstream>
#include <map>
#include <set>
#include <limits>
#include <unordered_set>
#include <vector>
#include <utility>
#include <stdexcept>
#include <string>
#include <iostream>
#include <cstdlib>
#include <iomanip>
#include <nlohmann/json.hpp>

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


    long long ParseGmshInteger( const std::string & token, const std::string & field )
    {
        std::size_t consumed = 0;
        long long value = 0;
        try
        {
            value = std::stoll( token, &consumed );
        }
        catch ( const std::exception & )
        {
            throw std::runtime_error( "Invalid Gmsh " + field + ": " + token );
        }
        if ( consumed != token.size() )
            throw std::runtime_error( "Invalid Gmsh " + field + ": " + token );
        return value;
    }

    double ParseGmshReal( const std::string & token, const std::string & field )
    {
        std::size_t consumed = 0;
        double value = 0.0;
        try
        {
            value = std::stod( token, &consumed );
        }
        catch ( const std::exception & )
        {
            throw std::runtime_error( "Invalid Gmsh " + field + ": " + token );
        }
        if ( consumed != token.size() )
            throw std::runtime_error( "Invalid Gmsh " + field + ": " + token );
        return value;
    }

    int GmshElementNodeCount( int type )
    {
        // Gmsh MSH element types 1-31: common line, surface and volume
        // elements through second order, including incomplete variants.
        switch ( type )
        {
            case 1: return 2;
            case 2: return 3;
            case 3: return 4;
            case 4: return 4;
            case 5: return 8;
            case 6: return 6;
            case 7: return 5;
            case 8: return 3;
            case 9: return 6;
            case 10: return 9;
            case 11: return 10;
            case 12: return 27;
            case 13: return 18;
            case 14: return 14;
            case 15: return 1;
            case 16: return 8;
            case 17: return 20;
            case 18: return 15;
            case 19: return 13;
            case 20: return 9;
            case 21: return 10;
            case 22: return 12;
            case 23: return 15;
            case 24: return 15;
            case 25: return 21;
            case 26: return 4;
            case 27: return 5;
            case 28: return 6;
            case 29: return 20;
            case 30: return 35;
            case 31: return 56;
            default:
                throw std::runtime_error(
                    "Unsupported Gmsh element type " + std::to_string( type ) +
                    " in MSH 4.1 ASCII reader" );
        }
    }

    struct GmshMeshSummary
    {
        long long nodeCount{ 0 };
        long long elementCount{ 0 };
        int versionMajor{ 0 };
        int versionMinor{ 0 };
        std::string versionText;
        std::string elementTypeCounts;
    };

    GmshMeshSummary ReadGmsh41Ascii( const std::string & content )
    {
        std::istringstream stream( content );
        std::vector< std::string > tokens;
        std::string token;
        while ( stream >> token ) tokens.push_back( token );

        auto findSection = [&tokens]( const std::string & name ) -> std::size_t
        {
            for ( std::size_t i = 0; i < tokens.size(); ++i )
                if ( tokens[ i ] == name ) return i + 1;
            throw std::runtime_error( "Missing Gmsh section " + name );
        };
        auto nextInteger = [&tokens]( std::size_t & pos, const std::string & field ) -> long long
        {
            if ( pos >= tokens.size() ) throw std::runtime_error( "Truncated Gmsh " + field );
            return ParseGmshInteger( tokens[ pos++ ], field );
        };
        auto nextReal = [&tokens]( std::size_t & pos, const std::string & field ) -> double
        {
            if ( pos >= tokens.size() ) throw std::runtime_error( "Truncated Gmsh " + field );
            return ParseGmshReal( tokens[ pos++ ], field );
        };

        GmshMeshSummary summary;
        summary.versionMajor = 4;
        summary.versionMinor = 1;
        summary.versionText = "4.1";

        std::unordered_set< long long > nodeTags;
        std::size_t pos = findSection( "$Nodes" );
        const long long nodeBlocks = nextInteger( pos, "node block count" );
        summary.nodeCount = nextInteger( pos, "node count" );
        const long long minNodeTag = nextInteger( pos, "minimum node tag" );
        const long long maxNodeTag = nextInteger( pos, "maximum node tag" );
        if ( nodeBlocks < 0 || summary.nodeCount < 1 ||
             minNodeTag < 1 || maxNodeTag < minNodeTag )
            throw std::runtime_error( "Invalid Gmsh 4.1 node section header" );

        long long observedNodes = 0;
        long long observedMinNodeTag = std::numeric_limits< long long >::max();
        long long observedMaxNodeTag = 0;
        for ( long long block = 0; block < nodeBlocks; ++block )
        {
            const long long entityDim = nextInteger( pos, "node entity dimension" );
            const long long entityTag = nextInteger( pos, "node entity tag" );
            const long long parametric = nextInteger( pos, "node parametric flag" );
            const long long blockNodes = nextInteger( pos, "node block size" );
            if ( entityDim < 0 || entityDim > 3 || entityTag < 1 ||
                 ( parametric != 0 && parametric != 1 ) || blockNodes < 0 )
                throw std::runtime_error( "Invalid Gmsh 4.1 node block header" );

            for ( long long i = 0; i < blockNodes; ++i )
            {
                const long long tag = nextInteger( pos, "node tag" );
                if ( tag < 1 || ! nodeTags.insert( tag ).second )
                    throw std::runtime_error( "Invalid or duplicate Gmsh node tag" );
                observedMinNodeTag = std::min( observedMinNodeTag, tag );
                observedMaxNodeTag = std::max( observedMaxNodeTag, tag );
            }
            for ( long long i = 0; i < blockNodes; ++i )
            {
                for ( int axis = 0; axis < 3; ++axis )
                {
                    if ( ! std::isfinite( nextReal( pos, "node coordinate" ) ) )
                        throw std::runtime_error( "Non-finite Gmsh node coordinate" );
                }
                if ( parametric )
                {
                    for ( long long axis = 0; axis < entityDim; ++axis )
                    {
                        if ( ! std::isfinite( nextReal( pos, "parametric node coordinate" ) ) )
                            throw std::runtime_error( "Non-finite Gmsh parametric coordinate" );
                    }
                }
            }
            observedNodes += blockNodes;
        }
        if ( observedNodes != summary.nodeCount ||
             observedMinNodeTag != minNodeTag || observedMaxNodeTag != maxNodeTag )
            throw std::runtime_error( "Gmsh 4.1 node block total does not match node count" );
        if ( pos >= tokens.size() || tokens[ pos ] != "$EndNodes" )
            throw std::runtime_error( "Malformed Gmsh 4.1 $Nodes section terminator" );

        pos = findSection( "$Elements" );
        const long long elementBlocks = nextInteger( pos, "element block count" );
        summary.elementCount = nextInteger( pos, "element count" );
        const long long minElementTag = nextInteger( pos, "minimum element tag" );
        const long long maxElementTag = nextInteger( pos, "maximum element tag" );
        if ( elementBlocks < 0 || summary.elementCount < 1 ||
             minElementTag < 1 || maxElementTag < minElementTag )
            throw std::runtime_error( "Invalid Gmsh 4.1 element section header" );

        std::map< int, long long > elementTypes;
        std::unordered_set< long long > elementTags;
        long long observedElements = 0;
        long long observedMinElementTag = std::numeric_limits< long long >::max();
        long long observedMaxElementTag = 0;
        for ( long long block = 0; block < elementBlocks; ++block )
        {
            const long long entityDim = nextInteger( pos, "element entity dimension" );
            const long long entityTag = nextInteger( pos, "element entity tag" );
            const int type = static_cast< int >( nextInteger( pos, "element type" ) );
            const long long blockElements = nextInteger( pos, "element block size" );
            if ( entityDim < 0 || entityDim > 3 || entityTag < 1 || blockElements < 0 )
                throw std::runtime_error( "Invalid Gmsh 4.1 element block header" );
            const int nodesPerElement = GmshElementNodeCount( type );
            for ( long long i = 0; i < blockElements; ++i )
            {
                const long long elementTag = nextInteger( pos, "element tag" );
                if ( elementTag < 1 || ! elementTags.insert( elementTag ).second )
                    throw std::runtime_error( "Invalid or duplicate Gmsh element tag" );
                observedMinElementTag = std::min( observedMinElementTag, elementTag );
                observedMaxElementTag = std::max( observedMaxElementTag, elementTag );
                for ( int node = 0; node < nodesPerElement; ++node )
                {
                    const long long nodeTag = nextInteger( pos, "element node tag" );
                    if ( nodeTags.find( nodeTag ) == nodeTags.end() )
                        throw std::runtime_error( "Gmsh element references an unknown node tag" );
                }
                ++observedElements;
                ++elementTypes[ type ];
            }
        }
        if ( observedElements != summary.elementCount ||
             observedMinElementTag != minElementTag || observedMaxElementTag != maxElementTag )
            throw std::runtime_error( "Gmsh 4.1 element block totals or tag range do not match the section header" );
        if ( pos >= tokens.size() || tokens[ pos ] != "$EndElements" )
            throw std::runtime_error( "Malformed Gmsh 4.1 $Elements section terminator" );

        std::ostringstream types;
        bool first = true;
        for ( const auto & entry : elementTypes )
        {
            if ( ! first ) types << ",";
            types << entry.first << ":" << entry.second;
            first = false;
        }
        summary.elementTypeCounts = types.str();
        return summary;
    }


    int GmshToVtkCellType( int gmshType )
    {
        switch ( gmshType )
        {
            case 1: return VTK_TYPE::LINE;
            case 2: return VTK_TYPE::TRIANGLE;
            case 3: return VTK_TYPE::QUADRILATERAL;
            case 4: return VTK_TYPE::TETRAHEDRON;
            case 5: return VTK_TYPE::HEXAHEDRON;
            case 6: return VTK_TYPE::PRISM;
            case 7: return VTK_TYPE::PYRAMID;
            default:
                throw std::runtime_error(
                    "Gmsh-to-CGNS conversion currently supports first-order "
                    "LINE_2, TRI_3, QUAD_4, TETRA_4, HEXA_8, PENTA_6 and PYRA_5 cells; "
                    "unsupported Gmsh element type " + std::to_string( gmshType ) );
        }
    }

    int GmshPhysicalNameToCgnsBcType( const std::string & physicalName )
    {
        std::string key;
        for ( unsigned char ch : physicalName )
        {
            if ( std::isalnum( ch ) )
                key.push_back( static_cast< char >( std::tolower( ch ) ) );
        }

        if ( key == "inflow" )
            return static_cast< int >( CGNS_ENUMV( BCInflow ) );
        if ( key == "outflow" )
            return static_cast< int >( CGNS_ENUMV( BCOutflow ) );
        if ( key == "symmetry" || key == "symmetryplane" )
            return static_cast< int >( CGNS_ENUMV( BCSymmetryPlane ) );
        if ( key == "farfield" )
            return static_cast< int >( CGNS_ENUMV( BCFarfield ) );
        if ( key == "wall" || key == "viscouswall" )
            return static_cast< int >( CGNS_ENUMV( BCWall ) );

        return static_cast< int >( CGNS_ENUMV( BCTypeUserDefined ) );
    }

    const char * CgnsBcTypeName( int bcType )
    {
        if ( bcType == static_cast< int >( CGNS_ENUMV( BCInflow ) ) ) return "BCInflow";
        if ( bcType == static_cast< int >( CGNS_ENUMV( BCOutflow ) ) ) return "BCOutflow";
        if ( bcType == static_cast< int >( CGNS_ENUMV( BCSymmetryPlane ) ) ) return "BCSymmetryPlane";
        if ( bcType == static_cast< int >( CGNS_ENUMV( BCFarfield ) ) ) return "BCFarfield";
        if ( bcType == static_cast< int >( CGNS_ENUMV( BCWall ) ) ) return "BCWall";
        return "BCTypeUserDefined";
    }

    UnstructuredMeshData BuildUnstructuredMeshFromGmsh( const GmshMesh & mesh, int & cellDimension )
    {
        if ( mesh.nodes.empty() || mesh.elements.empty() )
            throw std::runtime_error( "Gmsh mesh has no nodes or elements" );
        cellDimension = 0;
        for ( const GmshElement & element : mesh.elements )
            cellDimension = std::max( cellDimension, element.entityDimension );
        if ( cellDimension < 1 || cellDimension > 3 )
            throw std::runtime_error( "Gmsh mesh has no supported cell dimension" );
        if ( mesh.nodes.size() > static_cast< size_t >( std::numeric_limits< int >::max() ) )
            throw std::runtime_error( "Gmsh node count exceeds OneFLOW index range" );

        UnstructuredMeshData grid;
        grid.ndim = cellDimension;
        grid.nPoin = static_cast< int >( mesh.nodes.size() );
        grid.xN.reserve( mesh.nodes.size() );
        grid.yN.reserve( mesh.nodes.size() );
        grid.zN.reserve( mesh.nodes.size() );
        for ( const GmshNode & node : mesh.nodes )
        {
            grid.xN.push_back( node.coordinates[ 0 ] );
            grid.yN.push_back( node.coordinates[ 1 ] );
            grid.zN.push_back( node.coordinates[ 2 ] );
        }
        std::map< std::pair<int, long long>, std::vector< std::pair<int, IntField> > > boundaryGroups;
        for ( const GmshElement & element : mesh.elements )
        {
            if ( element.entityDimension != cellDimension &&
                 element.entityDimension != cellDimension - 1 ) continue;
            const int vtkType = GmshToVtkCellType( element.type );
            IntField connectivity;
            connectivity.reserve( element.nodeTags.size() );
            for ( long long nodeTag : element.nodeTags )
            {
                const auto iter = mesh.nodeTagToIndex.find( nodeTag );
                if ( iter == mesh.nodeTagToIndex.end() )
                    throw std::runtime_error( "Gmsh element references an unmapped source node tag" );
                if ( iter->second > static_cast< size_t >( std::numeric_limits< int >::max() ) )
                    throw std::runtime_error( "Gmsh node index exceeds OneFLOW index range" );
                connectivity.push_back( static_cast< int >( iter->second ) );
            }
            if ( element.entityDimension == cellDimension )
            {
                grid.elems.push_back( std::move( connectivity ) );
                grid.elemVTKType.push_back( vtkType );
                continue;
            }
            const auto entity = mesh.entityPhysicalTags.find(
                std::make_pair(element.entityDimension, element.entityTag));
            if (entity == mesh.entityPhysicalTags.end()) continue;
            for (long long physicalTag : entity->second)
            {
                const auto groupKey = std::make_pair(element.entityDimension, physicalTag);
                boundaryGroups[groupKey].push_back(std::make_pair(vtkType, connectivity));
            }
        }
        grid.mmark.CreateMarkerList(static_cast<int>(boundaryGroups.size()));
        grid.mmark.nMarker = static_cast<int>(boundaryGroups.size());
        int markerIndex = 0;
        for (const auto & group : boundaryGroups)
        {
            Marker & marker = grid.mmark.markerList[markerIndex++];
            const auto name = mesh.physicalNames.find(group.first);
            marker.name = name == mesh.physicalNames.end()
                ? "PhysicalGroup_" + std::to_string(group.first.first) + "_" +
                  std::to_string(group.first.second)
                : name->second;
            if (marker.name.size() > 32)
                throw std::runtime_error("Gmsh physical group name exceeds the CGNS 32-character limit: " + marker.name);
            marker.bcName = "GmshPhysicalGroup";
            marker.cgns_bcType = GmshPhysicalNameToCgnsBcType( marker.name );
            for (const auto & boundaryElement : group.second)
            {
                marker.eTypes.push_back(boundaryElement.first);
                marker.elems.push_back(boundaryElement.second);
            }
            marker.nElem = static_cast<int>(marker.elems.size());
        }
        if ( grid.elems.empty() )
            throw std::runtime_error( "Gmsh mesh contains no supported top-dimensional cells" );
        if ( grid.elems.size() > static_cast< size_t >( std::numeric_limits< int >::max() ) )
            throw std::runtime_error( "Gmsh cell count exceeds OneFLOW index range" );
        grid.nElem = static_cast< int >( grid.elems.size() );
        return grid;
    }

    void WriteGmshMeshAsCgns(
        const GmshMesh & mesh,
        int & cellDimension,
        const std::filesystem::path & target )
    {
        UnstructuredMeshData grid = BuildUnstructuredMeshFromGmsh( mesh, cellDimension );
        CgnsFactory factory;
        CgnsZbase & cgnsZbase = *factory.cgnsZbase;
        cgnsZbase.nBases = 1;
        cgnsZbase.InitCgnsBase();
        CgnsBase * cgnsBase = cgnsZbase.GetCgnsBase( 0 );
        cgnsBase->baseName = "Base";
        cgnsBase->celldim = cellDimension;
        cgnsBase->phydim = 3;
        cgnsBase->nZones = 1;
        cgnsBase->AllocateAllCgnsZones();
        CgnsZone * cgnsZone = cgnsBase->GetCgnsZone( 0 );
        UnstructuredCgnsZoneWriter::Write(grid, *cgnsZone, false);
        factory.WriteCgnsFile( target.string() );
    }

    void ConvertGmshToCgns(
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

        std::ifstream input( source, std::ios::binary );
        if ( ! input ) throw std::runtime_error( "Cannot open Gmsh input: " + source.string() );
        const std::string content{ std::istreambuf_iterator< char >( input ),
                                   std::istreambuf_iterator< char >() };
        std::istringstream lines( content );
        std::string line;
        GmshMeshSummary summary;
        GmshMesh parsedMesh;
        bool hasMaterializedMesh = false;
        bool hasFormat = false;
        while ( std::getline( lines, line ) )
        {
            if ( ! line.empty() && line.back() == '\r' ) line.pop_back();
            if ( line == "$MeshFormat" )
            {
                if ( ! std::getline( lines, line ) )
                    throw std::runtime_error( "Truncated Gmsh format header" );
                if ( ! line.empty() && line.back() == '\r' ) line.pop_back();
                std::istringstream row( line );
                double version = 0.0;
                int fileType = -1, dataSize = 0;
                if ( ! ( row >> version >> fileType >> dataSize ) || fileType != 0 )
                    throw std::runtime_error( "Gmsh reader supports ASCII MSH only" );
                if ( version >= 4.1 && version < 4.2 )
                {
                    summary = ReadGmsh41Ascii( content );
                    parsedMesh = ReadGmsh41AsciiMesh( content );
                    if ( static_cast<long long>(parsedMesh.nodes.size()) != summary.nodeCount ||
                         static_cast<long long>(parsedMesh.elements.size()) != summary.elementCount )
                        throw std::runtime_error( "Gmsh value reader count mismatch" );
                    hasMaterializedMesh = true;
                    (void)BuildGmshMeshReadReport( parsedMesh );
                    if ( ! target.parent_path().empty() )
                        std::filesystem::create_directories( target.parent_path() );
                }
                else if ( version >= 2.0 && version < 3.0 )
                {
                    summary.versionMajor = 2;
                    summary.versionMinor = static_cast< int >( ( version - 2.0 ) * 10.0 + 0.5 );
                    summary.versionText = line.substr( 0, line.find_first_of( " \t" ) );
                    // Keep the established MSH 2.x line parser below.
                    std::istringstream oldLines( content );
                    while ( std::getline( oldLines, line ) )
                    {
                        if ( ! line.empty() && line.back() == '\r' ) line.pop_back();
                        if ( line == "$Nodes" )
                        {
                            if ( ! std::getline( oldLines, line ) ) throw std::runtime_error( "Truncated Gmsh nodes section" );
                            std::istringstream count( line );
                            if ( ! ( count >> summary.nodeCount ) || summary.nodeCount < 1 ) throw std::runtime_error( "Invalid Gmsh node count" );
                            for ( long long i = 0; i < summary.nodeCount; ++i )
                            {
                                if ( ! std::getline( oldLines, line ) ) throw std::runtime_error( "Truncated Gmsh node records" );
                                std::istringstream node( line );
                                long long id = 0; double x = 0.0, y = 0.0, z = 0.0;
                                if ( ! ( node >> id >> x >> y >> z ) || id < 1 ) throw std::runtime_error( "Invalid Gmsh node record" );
                            }
                        }
                        else if ( line == "$Elements" )
                        {
                            if ( ! std::getline( oldLines, line ) ) throw std::runtime_error( "Truncated Gmsh elements section" );
                            std::istringstream count( line );
                            if ( ! ( count >> summary.elementCount ) || summary.elementCount < 1 ) throw std::runtime_error( "Invalid Gmsh element count" );
                            for ( long long i = 0; i < summary.elementCount; ++i )
                            {
                                if ( ! std::getline( oldLines, line ) ) throw std::runtime_error( "Truncated Gmsh element records" );
                                std::istringstream element( line );
                                long long id = 0, ignored = 0; int type = 0, tags = 0;
                                if ( ! ( element >> id >> type >> tags ) || tags < 0 ) throw std::runtime_error( "Invalid Gmsh element record" );
                                for ( int tag = 0; tag < tags; ++tag )
                                    if ( ! ( element >> ignored ) ) throw std::runtime_error( "Truncated Gmsh element tags" );
                                int nodes = 0;
                                while ( element >> ignored ) ++nodes;
                                if ( nodes == 0 ) throw std::runtime_error( "Gmsh element has no node references" );
                            }
                        }
                    }
                }
                else
                {
                    throw std::runtime_error( "Supported Gmsh ASCII versions are MSH 4.1 and MSH 2.x" );
                }
                hasFormat = true;
                break;
            }
        }
        if ( ! hasFormat || summary.nodeCount < 1 || summary.elementCount < 1 )
            throw std::runtime_error( "Gmsh input requires valid $MeshFormat, $Nodes and $Elements sections" );

        if ( ! target.parent_path().empty() )
            std::filesystem::create_directories( target.parent_path() );
        if ( ! hasMaterializedMesh )
            throw std::runtime_error(
                "MSH 2.x input is validated but not yet materialized for CGNS conversion" );
        int cellDimension = 0;
        WriteGmshMeshAsCgns( parsedMesh, cellDimension, target );

        // Reports belong to the case result area, not beside source files or
        // in the repository. Keep this output local to the conversion case.
        const std::filesystem::path resultDirectory =
            std::filesystem::path( caseDir ) / "results";
        std::filesystem::create_directories( resultDirectory );
        const std::string reportStem = target.stem().string();
        const std::filesystem::path jsonReportPath =
            resultDirectory / ( reportStem + ".gmsh-read-report.json" );
        std::ofstream jsonReport( jsonReportPath, std::ios::trunc );
        if ( ! jsonReport )
            throw std::runtime_error( "Cannot create Gmsh JSON report: " + jsonReportPath.string() );
        jsonReport << GmshMeshReadReportToJson( BuildGmshMeshReadReport( parsedMesh ) );
        if ( ! jsonReport )
            throw std::runtime_error( "Failed writing Gmsh JSON report: " + jsonReportPath.string() );

        const std::filesystem::path markdownReportPath =
            resultDirectory / ( reportStem + ".gmsh-report.md" );
        std::ofstream markdownReport( markdownReportPath, std::ios::trunc );
        if ( ! markdownReport )
            throw std::runtime_error( "Cannot create Gmsh Markdown report: " + markdownReportPath.string() );
        markdownReport << GmshMeshReportToMarkdown( parsedMesh, target.filename().string() );
        markdownReport << "\n## CGNS Boundary Condition Mapping\n\n"
                       << "| Physical Group | CGNS BCType_t |\n|---|---|\n";
        std::map< std::pair<int, long long>, std::string > boundaryPhysicalGroups;
        for ( const GmshElement & element : parsedMesh.elements )
        {
            if ( element.entityDimension != cellDimension - 1 ) continue;
            const auto entity = parsedMesh.entityPhysicalTags.find(
                std::make_pair( element.entityDimension, element.entityTag ) );
            if ( entity == parsedMesh.entityPhysicalTags.end() ) continue;
            for ( long long physicalTag : entity->second )
            {
                const auto groupKey = std::make_pair( element.entityDimension, physicalTag );
                const auto name = parsedMesh.physicalNames.find( groupKey );
                boundaryPhysicalGroups[groupKey] = name == parsedMesh.physicalNames.end()
                    ? "PhysicalGroup_" + std::to_string( groupKey.first ) + "_" +
                      std::to_string( groupKey.second )
                    : name->second;
            }
        }
        for ( const auto & group : boundaryPhysicalGroups )
        {
            const int bcType = GmshPhysicalNameToCgnsBcType( group.second );
            markdownReport << "| " << group.second << " | "
                           << CgnsBcTypeName( bcType ) << " |\n";
        }
        if ( ! markdownReport )
            throw std::runtime_error( "Failed writing Gmsh Markdown report: " + markdownReportPath.string() );

        std::cout << "Gmsh MSH 4.1 conversion completed: "
                  << parsedMesh.nodes.size() << " nodes, "
                  << parsedMesh.elements.size() << " elements." << std::endl;
        std::cout << "Gmsh reports: " << jsonReportPath.string()
                  << " and " << markdownReportPath.string() << std::endl;
    }

    void ConvertGmshToOneFLOW(
        const GridConfig & config,
        const std::string & caseDir )
    {
        // Keep CGNS as the explicit, independently readable interchange
        // artifact. This route composes the existing Gmsh writer and CGNS
        // reader instead of maintaining a second Gmsh-to-runtime mapping.
        const std::filesystem::path sourceName( config.sourceFile );
        const std::string interchangeName =
            sourceName.stem().string() + ".gmsh-to-oneflow.cgns";

        GridConfig cgnsConfig = config;
        cgnsConfig.targetType = GridFileType::CGNS;
        cgnsConfig.targetFile =
            ( std::filesystem::path( "results" ) / interchangeName ).generic_string();
        ConvertGmshToCgns( cgnsConfig, caseDir );

        GridConfig oneFlowConfig = cgnsConfig;
        oneFlowConfig.sourceType = GridFileType::CGNS;
        oneFlowConfig.sourceFile = cgnsConfig.targetFile;
        oneFlowConfig.targetType = GridFileType::OneFLOW;
        oneFlowConfig.topology = GridTopology::Unstructured;

        // The CGNS interchange file is an input to the OneFLOW reader, not
        // the native OneFLOW output. Never let CalcGrid::Dump overwrite it.
        // Respect an explicitly configured native-grid target; otherwise
        // keep the generated OneFLOW grid beside the interchange evidence.
        oneFlowConfig.targetFile = config.targetFile;
        if ( oneFlowConfig.targetFile.empty() )
        {
            const std::string nativeGridName =
                sourceName.stem().string() + ".oneflow.grid";
            oneFlowConfig.targetFile =
                ( std::filesystem::absolute( caseDir ) /
                  "results" / nativeGridName ).lexically_normal().string();
        }

        CgnsFactory cgnsFactory;
        cgnsFactory.GenerateGrid( oneFlowConfig, caseDir );
    }


    UnstructuredMeshData BuildUnstructuredMeshFromStl(
        const StlMesh & mesh,
        Real weldTolerance )
    {
        if ( mesh.triangles.empty() )
            throw std::runtime_error( "STL input contains no triangles" );
        if ( ! std::isfinite( weldTolerance ) || weldTolerance <= 0.0 )
            throw std::invalid_argument( "stlWeldTolerance must be finite and greater than zero" );
        UnstructuredMeshData grid;
        grid.ndim = 2;

        // The ADT supplies candidates only; PointLocator applies the exact
        // Euclidean tolerance test before reusing a node.
        RealField pmin( 3 ), pmax( 3 );
        for ( size_t axis = 0; axis < 3; ++axis )
        {
            pmin[ axis ] = std::numeric_limits< Real >::max();
            pmax[ axis ] = -std::numeric_limits< Real >::max();
        }
        for ( const StlTriangle & triangle : mesh.triangles )
        {
            for ( const auto & vertex : triangle.vertices )
            {
                for ( size_t axis = 0; axis < 3; ++axis )
                {
                    pmin[ axis ] = std::min( pmin[ axis ], static_cast< Real >( vertex[ axis ] ) );
                    pmax[ axis ] = std::max( pmax[ axis ], static_cast< Real >( vertex[ axis ] ) );
                }
            }
        }

        PointLocator vertexLocator;
        vertexLocator.Initialize( pmin, pmax, static_cast< Real >( weldTolerance ) );
        for ( const StlTriangle & triangle : mesh.triangles )
        {
            IntField connectivity;
            connectivity.reserve( 3 );
            for ( const auto & vertex : triangle.vertices )
            {
                const int pointId = vertexLocator.AddPoint(
                    static_cast< Real >( vertex[ 0 ] ),
                    static_cast< Real >( vertex[ 1 ] ),
                    static_cast< Real >( vertex[ 2 ] ) );
                if ( pointId < 0 )
                    throw std::runtime_error( "STL vertex locator returned an invalid point ID" );
                if ( static_cast< size_t >( pointId ) == grid.xN.size() )
                {
                    if ( grid.xN.size() >= static_cast< size_t >( std::numeric_limits< int >::max() ) )
                        throw std::runtime_error( "STL vertex count exceeds OneFLOW index range" );
                    grid.xN.push_back( vertex[ 0 ] );
                    grid.yN.push_back( vertex[ 1 ] );
                    grid.zN.push_back( vertex[ 2 ] );
                }
                connectivity.push_back( pointId );
            }
            if ( connectivity[ 0 ] == connectivity[ 1 ] ||
                 connectivity[ 1 ] == connectivity[ 2 ] ||
                 connectivity[ 2 ] == connectivity[ 0 ] )
                throw std::runtime_error( "STL triangle collapses under stlWeldTolerance" );
            grid.elems.push_back( std::move( connectivity ) );
            grid.elemVTKType.push_back( VTK_TYPE::TRIANGLE );
        }
        grid.nPoin = static_cast< int >( grid.xN.size() );
        grid.nElem = static_cast< int >( grid.elems.size() );
        grid.mmark.CreateMarkerList( 0 );
        grid.mmark.nMarker = 0;
        return grid;
    }


    void ConvertStlExternalFlowToGmsh(
        const GridConfig& config,
        const std::string& caseDir,
        const nlohmann::json& options,
        const std::filesystem::path& source,
        const std::filesystem::path& target)
    {
        using Json = nlohmann::json;
        namespace fs = std::filesystem;

        const std::string surfaceName = options.value("surface_name", std::string("fluid"));
        const std::string bodyName = options.value("body_boundary_name", std::string("wall"));
        const double bodyMeshSize = options.value("mesh_size", 0.15);
        const double outerMeshSize = options.value("outer_mesh_size", bodyMeshSize * 2.5);
        if (surfaceName.empty() || bodyName.empty() ||
            surfaceName.size() > 32 || bodyName.size() > 32)
            throw std::invalid_argument("External-flow physical names must contain 1 to 32 characters");
        const auto validName = [](const std::string& name)
        {
            return !name.empty() && std::all_of(name.begin(), name.end(),
                [](unsigned char value)
                {
                    return std::isalnum(value) || value == '_' || value == '.' || value == '-';
                });
        };
        if (!validName(surfaceName) || !validName(bodyName))
            throw std::invalid_argument("External-flow physical names may contain only letters, digits, _, . and -");
        if (!std::isfinite(bodyMeshSize) || bodyMeshSize <= 0.0 ||
            !std::isfinite(outerMeshSize) || outerMeshSize <= 0.0)
            throw std::invalid_argument("External-flow mesh sizes must be positive and finite");

        const Json& domain = options.at("domain");
        const double xmin = domain.at("xmin").get<double>();
        const double xmax = domain.at("xmax").get<double>();
        const double ymin = domain.at("ymin").get<double>();
        const double ymax = domain.at("ymax").get<double>();
        if (!std::isfinite(xmin) || !std::isfinite(xmax) ||
            !std::isfinite(ymin) || !std::isfinite(ymax) ||
            xmax - xmin <= config.stlWeldTolerance ||
            ymax - ymin <= config.stlWeldTolerance)
            throw std::invalid_argument("External-flow domain bounds must be finite and have positive dimensions");

        const StlMesh stl = StlMeshReader::Read(source.string());
        const StlPlanarBoundary body =
            ExtractStlPlanarBoundary(stl, config.stlWeldTolerance);
        if (body.vertices.size() < 3)
            throw std::invalid_argument("External-flow STL body must have at least three boundary edges");
        for (const auto& point : body.vertices)
        {
            if (point[0] <= xmin + config.stlWeldTolerance ||
                point[0] >= xmax - config.stlWeldTolerance ||
                point[1] <= ymin + config.stlWeldTolerance ||
                point[1] >= ymax - config.stlWeldTolerance)
                throw std::invalid_argument("External-flow outer domain must strictly enclose the STL body");
        }

        const Json& boundaryOptions = domain.at("boundaries");
        const std::array<std::string, 4> sides{"left", "right", "bottom", "top"};
        std::map<std::string, std::vector<int>> physicalCurves;
        std::set<std::string> outerNames;
        for (const std::string& side : sides)
        {
            const std::string name = boundaryOptions.at(side).get<std::string>();
            if (!validName(name) || name.size() > 32 || name == bodyName)
                throw std::invalid_argument("Invalid or conflicting external-flow boundary name: " + name);
            outerNames.insert(name);
        }
        if (outerNames.size() < 3)
            throw std::invalid_argument("External-flow domain requires distinct inlet, outlet and farfield groups");

        const fs::path resultDirectory = fs::path(caseDir) / "results";
        fs::create_directories(resultDirectory);
        if (!target.parent_path().empty())
            fs::create_directories(target.parent_path());
        const fs::path geoPath =
            resultDirectory / (target.stem().string() + ".geo");
        const fs::path reportPath =
            resultDirectory / (target.stem().string() + ".stl-to-gmsh-report.json");

        std::ofstream geo(geoPath, std::ios::trunc);
        if (!geo)
            throw std::runtime_error("Cannot create Gmsh geometry: " + geoPath.string());
        geo << std::setprecision(std::numeric_limits<double>::max_digits10);
        geo << "SetFactory(\"Built-in\");\n"
            << "Mesh.MshFileVersion = 4.1;\n"
            << "Mesh.Binary = 0;\n"
            << "Mesh.ElementOrder = 1;\n"
            << "Mesh.RecombineAll = 0;\n"
            << "Mesh.Algorithm = 6;\n"
            << "Mesh.MeshSizeMin = " << std::min(bodyMeshSize, outerMeshSize) << ";\n"
            << "Mesh.MeshSizeMax = " << std::max(bodyMeshSize, outerMeshSize) << ";\n";

        const std::array<std::array<double, 3>, 4> outerPoints{{
            {xmin, ymin, 0.0}, {xmax, ymin, 0.0},
            {xmax, ymax, 0.0}, {xmin, ymax, 0.0}}};
        for (std::size_t i = 0; i < outerPoints.size(); ++i)
        {
            const auto& point = outerPoints[i];
            geo << "Point(" << i + 1 << ") = {"
                << point[0] << ", " << point[1] << ", " << point[2]
                << ", " << outerMeshSize << "};\n";
        }
        for (int i = 0; i < 4; ++i)
            geo << "Line(" << i + 1 << ") = {"
                << i + 1 << ", " << (i + 1) % 4 + 1 << "};\n";

        const std::array<std::pair<std::string, int>, 4> sideCurves{{
            {"bottom", 1}, {"right", 2}, {"top", 3}, {"left", 4}}};
        for (const auto& side : sideCurves)
            physicalCurves[boundaryOptions.at(side.first).get<std::string>()].push_back(side.second);

        const int firstBodyPoint = 5;
        const int firstBodyLine = 5;
        for (std::size_t i = 0; i < body.vertices.size(); ++i)
        {
            const auto& point = body.vertices[i];
            geo << "Point(" << firstBodyPoint + static_cast<int>(i) << ") = {"
                << point[0] << ", " << point[1] << ", " << point[2]
                << ", " << bodyMeshSize << "};\n";
        }
        std::vector<int> reversedBodyLoop;
        std::vector<int> bodyPhysicalLines;
        for (std::size_t i = 0; i < body.vertices.size(); ++i)
        {
            const int lineId = firstBodyLine + static_cast<int>(i);
            const int pointId = firstBodyPoint + static_cast<int>(i);
            const int nextPointId =
                firstBodyPoint + static_cast<int>((i + 1) % body.vertices.size());
            geo << "Line(" << lineId << ") = {"
                << pointId << ", " << nextPointId << "};\n";
            reversedBodyLoop.push_back(-lineId);
            bodyPhysicalLines.push_back(lineId);
        }
        std::reverse(reversedBodyLoop.begin(), reversedBodyLoop.end());
        geo << "Curve Loop(1) = {1, 2, 3, 4};\n"
            << "Curve Loop(2) = {";
        for (std::size_t i = 0; i < reversedBodyLoop.size(); ++i)
        {
            if (i != 0) geo << ", ";
            geo << reversedBodyLoop[i];
        }
        geo << "};\n"
            << "Plane Surface(1) = {1, 2};\n";

        for (const auto& group : physicalCurves)
        {
            geo << "Physical Curve(\"" << group.first << "\") = {";
            for (std::size_t i = 0; i < group.second.size(); ++i)
            {
                if (i != 0) geo << ", ";
                geo << group.second[i];
            }
            geo << "};\n";
        }
        geo << "Physical Curve(\"" << bodyName << "\") = {";
        for (std::size_t i = 0; i < bodyPhysicalLines.size(); ++i)
        {
            if (i != 0) geo << ", ";
            geo << bodyPhysicalLines[i];
        }
        geo << "};\n"
            << "Physical Surface(\"" << surfaceName << "\") = {1};\n";
        geo.flush();
        if (!geo)
            throw std::runtime_error("Failed writing Gmsh geometry: " + geoPath.string());
        geo.close();

        if (config.gmshExecutable.empty())
            throw std::invalid_argument("gmshExecutable must not be empty");
        const auto pathToUtf8 = [](const fs::path& path)
        {
            const auto value = path.u8string();
            return std::string(value.begin(), value.end());
        };
        const std::vector<std::string> gmshArguments{
            pathToUtf8(fs::absolute(geoPath)),
            "-2", "-format", "msh4", "-o",
            pathToUtf8(fs::absolute(target))};
        const ProcessResult process = RunProcess(
            config.gmshExecutable, gmshArguments,
            pathToUtf8(fs::absolute(fs::path(caseDir))));
        if (!process.started || process.exitCode != 0)
        {
            std::ostringstream diagnostic;
            diagnostic << "Gmsh external-flow mesh generation failed"
                << "\n  executable: " << config.gmshExecutable
                << "\n  working directory: " << pathToUtf8(fs::absolute(fs::path(caseDir)))
                << "\n  started: " << (process.started ? "true" : "false")
                << "\n  exit code: " << process.exitCode;
            if (!process.errorMessage.empty())
                diagnostic << "\n  process error: " << process.errorMessage;
            if (!process.output.empty())
                diagnostic << "\n  output:\n" << process.output;
            throw std::runtime_error(diagnostic.str());
        }

        std::ifstream meshInput(target, std::ios::binary);
        if (!meshInput)
            throw std::runtime_error("Gmsh did not create the requested mesh: " + target.string());
        const std::string meshContent{
            std::istreambuf_iterator<char>(meshInput),
            std::istreambuf_iterator<char>()};
        const GmshMesh generated = ReadGmsh41AsciiMesh(meshContent);
        const GmshMeshReadReport meshReport = BuildGmshMeshReadReport(generated);
        const double domainTolerance = std::max(
            config.stlWeldTolerance,
            64.0 * std::numeric_limits<double>::epsilon() *
                std::max({std::abs(xmin), std::abs(xmax),
                          std::abs(ymin), std::abs(ymax), 1.0}));
        if (std::abs(meshReport.minimumCoordinates[0] - xmin) > domainTolerance ||
            std::abs(meshReport.maximumCoordinates[0] - xmax) > domainTolerance ||
            std::abs(meshReport.minimumCoordinates[1] - ymin) > domainTolerance ||
            std::abs(meshReport.maximumCoordinates[1] - ymax) > domainTolerance ||
            std::abs(meshReport.maximumCoordinates[2] -
                     meshReport.minimumCoordinates[2]) > domainTolerance)
            throw std::runtime_error(
                "Generated external-flow mesh coordinates do not match the planar computational domain");

        std::size_t lineElements = 0;
        std::size_t triangleElements = 0;
        std::size_t quadrilateralElements = 0;
        for (const GmshElement& element : generated.elements)
        {
            if (element.type == 1) ++lineElements;
            else if (element.type == 2) ++triangleElements;
            else if (element.type == 3) ++quadrilateralElements;
        }
        if (generated.nodes.size() < 20 ||
            lineElements < body.vertices.size() + 4 ||
            triangleElements == 0 || quadrilateralElements != 0 ||
            meshReport.degenerateTriangleCount != 0)
            throw std::runtime_error(
                "External-flow mesh topology is incomplete: expected boundary lines and "
                "triangular fluid cells, observed " +
                std::to_string(generated.nodes.size()) + " nodes, " +
                std::to_string(lineElements) + " lines, " +
                std::to_string(triangleElements) + " triangles and " +
                std::to_string(quadrilateralElements) + " quadrilaterals");

        const Json qualityOptions = options.value("quality", Json::object());
        const double minimumQuality = qualityOptions.value(
            "minimum_triangle_quality", 0.01);
        const double minimumAngleDegrees = qualityOptions.value(
            "minimum_angle_degrees", 1.0);
        if (!std::isfinite(minimumQuality) || minimumQuality <= 0.0 ||
            minimumQuality > 1.0 || !std::isfinite(minimumAngleDegrees) ||
            minimumAngleDegrees <= 0.0 || minimumAngleDegrees >= 60.0)
            throw std::invalid_argument("External-flow triangle quality thresholds are invalid");

        double observedMinimumQuality = 1.0;
        double observedMinimumAngleDegrees = 60.0;
        for (const GmshElement& element : generated.elements)
        {
            if (element.type != 2 || element.nodeTags.size() != 3) continue;
            std::array<std::array<double, 3>, 3> points{};
            for (std::size_t i = 0; i < 3; ++i)
            {
                const auto found = generated.nodeTagToIndex.find(element.nodeTags[i]);
                if (found == generated.nodeTagToIndex.end())
                    throw std::runtime_error("TRI_3 references an unknown node while measuring quality");
                points[i] = generated.nodes[found->second].coordinates;
            }
            const auto distance = [](const auto& a, const auto& b)
            {
                const double x = a[0] - b[0];
                const double y = a[1] - b[1];
                const double z = a[2] - b[2];
                return std::sqrt(x * x + y * y + z * z);
            };
            const double a = distance(points[0], points[1]);
            const double b = distance(points[1], points[2]);
            const double c = distance(points[2], points[0]);
            const double twiceArea = 2.0 * (
                (points[1][0] - points[0][0]) * (points[2][1] - points[0][1]) -
                (points[1][1] - points[0][1]) * (points[2][0] - points[0][0]));
            const double area = 0.5 * std::abs(twiceArea);
            if (area <= 0.0 || a <= 0.0 || b <= 0.0 || c <= 0.0)
                throw std::runtime_error("External-flow mesh contains a degenerate TRI_3 element");
            const double quality = 4.0 * std::sqrt(3.0) * area /
                (a * a + b * b + c * c);
            const auto angle = [](double opposite, double adjacentA, double adjacentB)
            {
                const double cosine = (adjacentA * adjacentA + adjacentB * adjacentB -
                    opposite * opposite) / (2.0 * adjacentA * adjacentB);
                return std::acos(std::max(-1.0, std::min(1.0, cosine))) *
                    180.0 / 3.14159265358979323846;
            };
            observedMinimumQuality = std::min(observedMinimumQuality, quality);
            observedMinimumAngleDegrees = std::min({
                observedMinimumAngleDegrees,
                angle(b, a, c), angle(c, a, b), angle(a, b, c)});
        }
        if (observedMinimumQuality < minimumQuality ||
            observedMinimumAngleDegrees < minimumAngleDegrees)
            throw std::runtime_error(
                "External-flow mesh failed triangle quality thresholds: minimum quality " +
                std::to_string(observedMinimumQuality) + ", minimum angle " +
                std::to_string(observedMinimumAngleDegrees) + " degrees");

        std::set<std::pair<int, std::string>> expectedGroups;
        for (const auto& group : physicalCurves)
            expectedGroups.emplace(1, group.first);
        expectedGroups.emplace(1, bodyName);
        expectedGroups.emplace(2, surfaceName);
        std::set<std::pair<int, std::string>> actualGroups;
        for (const auto& physical : generated.physicalNames)
            actualGroups.emplace(physical.first.first, physical.second);
        if (actualGroups != expectedGroups)
            throw std::runtime_error("Generated external-flow MSH Physical Groups do not match the boundary contract");

        Json report;
        report["schema"] = "oneflow.stl-to-gmsh-validation-report";
        report["version"] = 1;
        report["case_type"] = "2d_external_flow";
        report["source_stl"] = source.lexically_normal().string();
        report["geometry"] = geoPath.lexically_normal().string();
        report["mesh"] = target.lexically_normal().string();
        report["gmsh_executable"] = config.gmshExecutable;
        report["gmsh_version"] = generated.version;
        report["stl_triangles"] = stl.triangles.size();
        report["stl_boundary_vertices"] = body.vertices.size();
        report["contract"] = {
            {"domain", {
                {"xmin", xmin}, {"xmax", xmax}, {"ymin", ymin}, {"ymax", ymax},
                {"boundaries", boundaryOptions}}},
            {"surface_name", surfaceName},
            {"body_boundary_name", bodyName},
            {"quality_limits", {
                {"minimum_triangle_quality", minimumQuality},
                {"minimum_angle_degrees", minimumAngleDegrees}}}};
        report["domain"] = report["contract"]["domain"];
        report["mesh_summary"] = {
            {"nodes", generated.nodes.size()},
            {"line_elements", lineElements},
            {"triangle_elements", triangleElements},
            {"quadrilateral_elements", quadrilateralElements},
            {"physical_group_count", meshReport.physicalGroupCount}};
        report["mesh_summary"]["triangle_quality"] = {
            {"minimum_quality", observedMinimumQuality},
            {"minimum_angle_degrees", observedMinimumAngleDegrees},
            {"minimum_area", meshReport.minimumTriangleArea},
            {"degenerate_count", meshReport.degenerateTriangleCount}};
        report["mesh_summary"]["physical_groups"] = Json::array();
        for (const auto& group : actualGroups)
            report["mesh_summary"]["physical_groups"].push_back({
                {"dimension", group.first}, {"name", group.second}});
        report["validated"] = true;
        std::ofstream reportOutput(reportPath, std::ios::trunc);
        if (!reportOutput)
            throw std::runtime_error("Cannot create STL-to-Gmsh report: " + reportPath.string());
        reportOutput << report.dump(2) << '\n';
        if (!reportOutput)
            throw std::runtime_error("Failed writing STL-to-Gmsh report: " + reportPath.string());

        std::cout << "STL external-flow mesh completed: "
            << generated.nodes.size() << " nodes, "
            << triangleElements << " TRI_3 cells, "
            << lineElements << " boundary LINE_2 elements." << std::endl;
    }

    void ConvertStlToGmsh(
        const GridConfig& config,
        const std::string& caseDir)
    {
        using Json = nlohmann::json;
        namespace fs = std::filesystem;

        if (config.layoutFile.empty())
            throw std::invalid_argument(
                "STL to Gmsh requires gridLayoutFileName to identify its JSON layout");

        const fs::path sourceName(config.sourceFile);
        const fs::path sourceRoot = config.ResolveSourceCaseDir(caseDir);
        const fs::path source = sourceName.is_absolute()
            ? sourceName : sourceRoot / sourceName;
        const fs::path layoutName(config.layoutFile);
        const fs::path layout = layoutName.is_absolute()
            ? layoutName : fs::path(caseDir) / layoutName;
        const fs::path targetName(config.targetFile);
        const fs::path target = targetName.is_absolute()
            ? targetName : fs::path(caseDir) / targetName;
        if (config.targetFile.empty())
            throw std::invalid_argument("STL to Gmsh requires targetGridFileName");

        std::ifstream layoutInput(layout, std::ios::binary);
        if (!layoutInput)
            throw std::runtime_error("Cannot open STL-to-Gmsh layout: " + layout.string());
        const Json options = Json::parse(layoutInput);

        if (options.contains("domain"))
        {
            ConvertStlExternalFlowToGmsh(config, caseDir, options, source, target);
            return;
        }

        const std::string surfaceName = options.value("surface_name", std::string("fluid_domain"));
        if (surfaceName != "fluid_domain")
            throw std::invalid_argument(
                "The initial STL-to-Gmsh contract requires surface_name = fluid_domain");
        const double meshSize = options.value("mesh_size", 1.0);
        if (!std::isfinite(meshSize) || meshSize <= 0.0)
            throw std::invalid_argument("STL-to-Gmsh mesh_size must be positive and finite");

        const Json cells = options.at("transfinite_cells");
        const int cellsX = cells.at("x").get<int>();
        const int cellsY = cells.at("y").get<int>();
        if (cellsX < 1 || cellsY < 1)
            throw std::invalid_argument("STL-to-Gmsh transfinite cell counts must be positive");

        const StlMesh stl = StlMeshReader::Read(source.string());
        const StlPlanarBoundary boundary =
            ExtractStlPlanarBoundary(stl, config.stlWeldTolerance);
        if (boundary.vertices.size() != 4)
            throw std::invalid_argument(
                "The first native STL-to-Gmsh stage requires a four-edge rectangular boundary");

        std::array<double, 4> bounds{
            boundary.vertices[0][0], boundary.vertices[0][0],
            boundary.vertices[0][1], boundary.vertices[0][1]};
        for (const auto& point : boundary.vertices)
        {
            bounds[0] = std::min(bounds[0], point[0]);
            bounds[1] = std::max(bounds[1], point[0]);
            bounds[2] = std::min(bounds[2], point[1]);
            bounds[3] = std::max(bounds[3], point[1]);
        }
        if (bounds[1] - bounds[0] <= config.stlWeldTolerance ||
            bounds[3] - bounds[2] <= config.stlWeldTolerance)
            throw std::invalid_argument("STL rectangular boundary has zero width or height");

        for (std::size_t i = 0; i < boundary.vertices.size(); ++i)
        {
            const auto& a = boundary.vertices[i];
            const auto& b = boundary.vertices[(i + 1) % boundary.vertices.size()];
            const bool horizontal = std::abs(a[1] - b[1]) <= config.stlWeldTolerance;
            const bool vertical = std::abs(a[0] - b[0]) <= config.stlWeldTolerance;
            if (horizontal == vertical)
                throw std::invalid_argument(
                    "STL boundary must be an axis-aligned rectangle for transfinite QUAD_4 generation");
        }

        const Json boundaryOptions = options.at("boundaries");
        if (!boundaryOptions.is_array() || boundaryOptions.size() != 4)
            throw std::invalid_argument(
                "STL-to-Gmsh layout must define exactly four boundary assignments");

        const std::array<std::string, 4> requiredNames{
            "inflow", "outflow", "symmetry", "farfield"};
        std::map<std::string, std::size_t> assignedCurves;
        std::set<std::size_t> assignedEdges;
        for (const Json& entry : boundaryOptions)
        {
            const std::string name = entry.at("name").get<std::string>();
            if (std::find(requiredNames.begin(), requiredNames.end(), name) ==
                requiredNames.end())
                throw std::invalid_argument("Unsupported outer boundary name: " + name);
            if (assignedCurves.count(name) != 0)
                throw std::invalid_argument("Duplicate STL-to-Gmsh boundary name: " + name);

            const Json endpoints = entry.at("endpoints");
            if (!endpoints.is_array() || endpoints.size() != 2)
                throw std::invalid_argument(name + " requires two [x, y] endpoints");
            std::array<std::size_t, 2> vertexIds{};
            for (std::size_t endpoint = 0; endpoint < 2; ++endpoint)
            {
                if (!endpoints[endpoint].is_array() || endpoints[endpoint].size() != 2)
                    throw std::invalid_argument(name + " endpoint must contain [x, y]");
                const double x = endpoints[endpoint][0].get<double>();
                const double y = endpoints[endpoint][1].get<double>();
                std::size_t matches = 0;
                for (std::size_t vertex = 0; vertex < boundary.vertices.size(); ++vertex)
                {
                    const auto& point = boundary.vertices[vertex];
                    if (std::hypot(point[0] - x, point[1] - y) <=
                        config.stlWeldTolerance)
                    {
                        vertexIds[endpoint] = vertex;
                        ++matches;
                    }
                }
                if (matches != 1)
                    throw std::invalid_argument(
                        name + " endpoint does not uniquely match an STL boundary vertex");
            }

            const std::size_t first = vertexIds[0];
            const std::size_t second = vertexIds[1];
            std::size_t edge = boundary.vertices.size();
            if ((first + 1) % boundary.vertices.size() == second)
                edge = first;
            else if ((second + 1) % boundary.vertices.size() == first)
                edge = second;
            if (edge == boundary.vertices.size() || !assignedEdges.insert(edge).second)
                throw std::invalid_argument(
                    name + " must identify one unique consecutive STL boundary edge");

            const auto& a = boundary.vertices[edge];
            const auto& b = boundary.vertices[(edge + 1) % boundary.vertices.size()];
            const bool isLeft =
                std::abs(a[0] - bounds[0]) <= config.stlWeldTolerance &&
                std::abs(b[0] - bounds[0]) <= config.stlWeldTolerance;
            const bool isRight =
                std::abs(a[0] - bounds[1]) <= config.stlWeldTolerance &&
                std::abs(b[0] - bounds[1]) <= config.stlWeldTolerance;
            const bool isBottom =
                std::abs(a[1] - bounds[2]) <= config.stlWeldTolerance &&
                std::abs(b[1] - bounds[2]) <= config.stlWeldTolerance;
            const bool isTop =
                std::abs(a[1] - bounds[3]) <= config.stlWeldTolerance &&
                std::abs(b[1] - bounds[3]) <= config.stlWeldTolerance;
            const bool matchesSide =
                (name == "inflow" && isLeft) ||
                (name == "outflow" && isRight) ||
                (name == "symmetry" && isBottom) ||
                (name == "farfield" && isTop);
            if (!matchesSide)
                throw std::invalid_argument(
                    name + " does not match its required left/right/bottom/top side");
            assignedCurves.emplace(name, edge + 1);
        }
        if (assignedEdges.size() != 4 || assignedCurves.size() != 4)
            throw std::invalid_argument("All four outer boundary edges must be assigned");

        const fs::path resultDirectory = fs::path(caseDir) / "results";
        fs::create_directories(resultDirectory);
        if (!target.parent_path().empty())
            fs::create_directories(target.parent_path());
        const fs::path geoPath =
            resultDirectory / (target.stem().string() + ".geo");
        const fs::path reportPath =
            resultDirectory / (target.stem().string() + ".stl-to-gmsh-report.json");

        std::ofstream geo(geoPath, std::ios::trunc);
        if (!geo)
            throw std::runtime_error("Cannot create Gmsh geometry: " + geoPath.string());
        geo << std::setprecision(std::numeric_limits<double>::max_digits10);
        geo << "SetFactory(\"Built-in\");\n"
            << "Mesh.MshFileVersion = 4.1;\n"
            << "Mesh.Binary = 0;\n"
            << "Mesh.ElementOrder = 1;\n"
            << "Mesh.RecombineAll = 1;\n";
        for (std::size_t i = 0; i < boundary.vertices.size(); ++i)
        {
            const auto& point = boundary.vertices[i];
            geo << "Point(" << i + 1 << ") = {"
                << point[0] << ", " << point[1] << ", " << point[2]
                << ", " << meshSize << "};\n";
        }
        for (std::size_t i = 0; i < boundary.vertices.size(); ++i)
            geo << "Line(" << i + 1 << ") = {"
                << i + 1 << ", " << (i + 1) % boundary.vertices.size() + 1 << "};\n";
        geo << "Curve Loop(1) = {1, 2, 3, 4};\n"
            << "Plane Surface(1) = {1};\n";
        for (std::size_t i = 0; i < boundary.vertices.size(); ++i)
        {
            const auto& a = boundary.vertices[i];
            const auto& b = boundary.vertices[(i + 1) % boundary.vertices.size()];
            const bool horizontal =
                std::abs(a[1] - b[1]) <= config.stlWeldTolerance;
            geo << "Transfinite Curve {" << i + 1 << "} = "
                << (horizontal ? cellsX + 1 : cellsY + 1) << ";\n";
        }
        geo << "Transfinite Surface {1};\n"
            << "Recombine Surface {1};\n";
        for (const auto& assignment : assignedCurves)
            geo << "Physical Curve(\"" << assignment.first << "\") = {"
                << assignment.second << "};\n";
        geo << "Physical Surface(\"" << surfaceName << "\") = {1};\n";
        geo.flush();
        if (!geo)
            throw std::runtime_error("Failed writing Gmsh geometry: " + geoPath.string());
        geo.close();

        if (config.gmshExecutable.empty())
            throw std::invalid_argument("gmshExecutable must not be empty");
        const auto pathToUtf8 = []( const fs::path & path )
        {
            const auto value = path.u8string();
            return std::string( value.begin(), value.end() );
        };
        const std::vector< std::string > gmshArguments{
            pathToUtf8( fs::absolute( geoPath ) ),
            "-2", "-format", "msh4", "-o",
            pathToUtf8( fs::absolute( target ) )
        };
        const ProcessResult process = RunProcess(
            config.gmshExecutable, gmshArguments,
            pathToUtf8( fs::absolute( fs::path( caseDir ) ) ) );
        if ( ! process.started || process.exitCode != 0 )
        {
            std::ostringstream diagnostic;
            diagnostic
                << "Gmsh process failed"
                << "\n  executable: " << config.gmshExecutable
                << "\n  working directory: "
                << pathToUtf8( fs::absolute( fs::path( caseDir ) ) )
                << "\n  started: " << ( process.started ? "true" : "false" )
                << "\n  exit code: " << process.exitCode;
            if ( ! process.errorMessage.empty() )
                diagnostic << "\n  process error: " << process.errorMessage;
            if ( ! process.output.empty() )
                diagnostic << "\n  output:\\n" << process.output;
            throw std::runtime_error( diagnostic.str() );
        }

        std::ifstream meshInput(target, std::ios::binary);
        if (!meshInput)
            throw std::runtime_error("Gmsh did not create the requested mesh: " + target.string());
        const std::string meshContent{
            std::istreambuf_iterator<char>(meshInput),
            std::istreambuf_iterator<char>()};
        const GmshMesh generated = ReadGmsh41AsciiMesh(meshContent);
        const GmshMeshReadReport meshReport = BuildGmshMeshReadReport(generated);

        std::size_t lineElements = 0;
        std::size_t quadrilateralElements = 0;
        std::size_t triangleElements = 0;
        for (const GmshElement& element : generated.elements)
        {
            if (element.type == 1) ++lineElements;
            else if (element.type == 2) ++triangleElements;
            else if (element.type == 3) ++quadrilateralElements;
        }
        const std::size_t expectedNodes =
            static_cast<std::size_t>(cellsX + 1) * static_cast<std::size_t>(cellsY + 1);
        const std::size_t expectedLines =
            2u * (static_cast<std::size_t>(cellsX) + static_cast<std::size_t>(cellsY));
        const std::size_t expectedQuads =
            static_cast<std::size_t>(cellsX) * static_cast<std::size_t>(cellsY);
        if (generated.nodes.size() != expectedNodes ||
            generated.elements.size() != expectedLines + expectedQuads ||
            lineElements != expectedLines ||
            quadrilateralElements != expectedQuads ||
            triangleElements != 0)
            throw std::runtime_error(
                "Generated MSH topology mismatch: expected " +
                std::to_string(expectedNodes) + " nodes, " +
                std::to_string(expectedLines) + " boundary lines, " +
                std::to_string(expectedQuads) + " QUAD_4 cells; observed " +
                std::to_string(generated.nodes.size()) + ", " +
                std::to_string(lineElements) + ", " +
                std::to_string(quadrilateralElements) + " quads and " +
                std::to_string(triangleElements) + " triangles");

        const std::set<std::pair<int, std::string>> expectedGroups{
            {1, "inflow"}, {1, "outflow"}, {1, "symmetry"},
            {1, "farfield"}, {2, surfaceName}};
        std::set<std::pair<int, std::string>> actualGroups;
        for (const auto& physical : generated.physicalNames)
            actualGroups.emplace(physical.first.first, physical.second);
        if (actualGroups != expectedGroups)
            throw std::runtime_error("Generated MSH Physical Groups do not match the boundary contract");

        Json report;
        report["schema"] = "oneflow.stl-to-gmsh-report";
        report["version"] = 1;
        report["source_stl"] = source.lexically_normal().string();
        report["geometry"] = geoPath.lexically_normal().string();
        report["mesh"] = target.lexically_normal().string();
        report["gmsh_executable"] = config.gmshExecutable;
        report["gmsh_version"] = generated.version;
        report["stl_triangles"] = stl.triangles.size();
        report["stl_boundary_vertices"] = boundary.vertices.size();
        report["transfinite_cells"] = {{"x", cellsX}, {"y", cellsY}};
        report["mesh_summary"] = {
            {"nodes", generated.nodes.size()},
            {"line_elements", lineElements},
            {"triangle_elements", triangleElements},
            {"quadrilateral_elements", quadrilateralElements},
            {"physical_group_count", meshReport.physicalGroupCount}};
        report["mesh_summary"]["physical_groups"] = Json::array();
        for (const auto& group : actualGroups)
        {
            report["mesh_summary"]["physical_groups"].push_back({
                {"dimension", group.first}, {"name", group.second}});
        }
        report["validated"] = true;
        std::ofstream reportOutput(reportPath, std::ios::trunc);
        if (!reportOutput)
            throw std::runtime_error("Cannot create STL-to-Gmsh report: " + reportPath.string());
        reportOutput << report.dump(2) << '\n';
        if (!reportOutput)
            throw std::runtime_error("Failed writing STL-to-Gmsh report: " + reportPath.string());

        std::cout << "STL to Gmsh completed: " << generated.nodes.size()
                  << " nodes, " << quadrilateralElements
                  << " QUAD_4 cells, " << lineElements
                  << " boundary LINE_2 elements." << std::endl;
    }

    void ConvertStlToCgns(
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
        const StlMesh mesh = StlMeshReader::Read( source.string() );
        UnstructuredMeshData grid = BuildUnstructuredMeshFromStl( mesh, config.stlWeldTolerance );
        if ( ! target.parent_path().empty() )
            std::filesystem::create_directories( target.parent_path() );
        CgnsFactory factory;
        CgnsZbase & cgnsZbase = *factory.cgnsZbase;
        cgnsZbase.nBases = 1;
        cgnsZbase.InitCgnsBase();
        CgnsBase * cgnsBase = cgnsZbase.GetCgnsBase( 0 );
        cgnsBase->baseName = "Base";
        cgnsBase->celldim = 2;
        cgnsBase->phydim = 3;
        cgnsBase->nZones = 1;
        cgnsBase->AllocateAllCgnsZones();
        CgnsZone * cgnsZone = cgnsBase->GetCgnsZone( 0 );
        UnstructuredCgnsZoneWriter::Write( grid, *cgnsZone, false );
        factory.WriteCgnsFile( target.string() );
        const std::filesystem::path resultDirectory = std::filesystem::path( caseDir ) / "results";
        std::filesystem::create_directories( resultDirectory );
        const std::filesystem::path reportPath =
            resultDirectory / ( target.stem().string() + ".stl-read-report.json" );
        std::ofstream report( reportPath, std::ios::trunc );
        if ( ! report ) throw std::runtime_error( "Cannot create STL report: " + reportPath.string() );
        report << "{\n"
               << "  \"format\": \"auto\",\n"
               << "  \"triangleCount\": " << mesh.triangles.size() << ",\n"
               << "  \"weldTolerance\": " << config.stlWeldTolerance << ",\n"
               << "  \"nodeCount\": " << grid.nPoin << ",\n"
               << "  \"cellCount\": " << grid.nElem << ",\n"
               << "  \"topology\": \"triangle_surface\"\n"
               << "}\n";
        if ( ! report ) throw std::runtime_error( "Failed writing STL report: " + reportPath.string() );
        std::cout << "STL to CGNS conversion completed: " << grid.nPoin
                  << " welded nodes, " << grid.nElem << " triangles." << std::endl;
    }

    void ConvertStl( const GridConfig & config, const std::string & caseDir )
    {
        switch ( config.targetType )
        {
            case GridFileType::CGNS:
                ConvertStlToCgns( config, caseDir );
                return;
            case GridFileType::Gmsh:
                ConvertStlToGmsh( config, caseDir );
                return;
            default:
                throw std::invalid_argument(
                    std::string( "Unsupported STL target type: " ) +
                    std::string( ToString( config.targetType ) ) );
        }
    }

    void ConvertGmsh( const GridConfig & config, const std::string & caseDir )
    {
        switch ( config.targetType )
        {
            case GridFileType::CGNS:
                ConvertGmshToCgns( config, caseDir );
                return;
            case GridFileType::OneFLOW:
                ConvertGmshToOneFLOW( config, caseDir );
                return;
            default:
                throw std::invalid_argument(
                    std::string( "Unsupported Gmsh target type: " ) +
                    std::string( ToString( config.targetType ) ) );
        }
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
        { GridFileType::STL,     &ConvertStl },
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
            return targetType == GridFileType::CGNS ||
                   targetType == GridFileType::OneFLOW;
        case GridFileType::STL:
            return targetType == GridFileType::CGNS ||
                   targetType == GridFileType::Gmsh;
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

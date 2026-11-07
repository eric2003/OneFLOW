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
#include "CgnsFactory.h"
#include "CgnsZbase.h"
#include "CgnsBase.h"
#include "CgnsZone.h"
#include "GridMediator.h"
#include "Su2Grid.h"
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
#include <limits>
#include <unordered_set>
#include <vector>
#include <utility>
#include <stdexcept>
#include <string>
#include <iostream>

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

        return static_cast< int >( CGNS_ENUMV( BCTypeUserDefined ) );
    }

    const char * CgnsBcTypeName( int bcType )
    {
        if ( bcType == static_cast< int >( CGNS_ENUMV( BCInflow ) ) ) return "BCInflow";
        if ( bcType == static_cast< int >( CGNS_ENUMV( BCOutflow ) ) ) return "BCOutflow";
        if ( bcType == static_cast< int >( CGNS_ENUMV( BCSymmetryPlane ) ) ) return "BCSymmetryPlane";
        if ( bcType == static_cast< int >( CGNS_ENUMV( BCFarfield ) ) ) return "BCFarfield";
        return "BCTypeUserDefined";
    }

    Su2Grid BuildSu2GridFromGmsh( const GmshMesh & mesh, int & cellDimension )
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

        Su2Grid grid;
        grid.ndim = cellDimension;
        grid.nPoin = static_cast< int >( mesh.nodes.size() );
        grid.nZone = 1;
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
                grid.elemId.push_back( static_cast< int >( grid.elemId.size() + 1 ) );
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
        Su2Grid grid = BuildSu2GridFromGmsh( mesh, cellDimension );
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
        grid.FillSU2CgnsZone( *cgnsZone, false );
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

    void ConvertGmsh( const GridConfig & config, const std::string & caseDir )
    {
        if ( config.targetType != GridFileType::CGNS )
            throw std::invalid_argument( "Initial Gmsh conversion supports targetGridType = cgns only" );
        ConvertGmshToCgns( config, caseDir );
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

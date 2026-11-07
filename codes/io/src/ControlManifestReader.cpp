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

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <set>
#include <string>
#include <system_error>
#include <vector>

namespace ONEFLOW {
namespace {

std::string Trim( const std::string & value )
{
    const std::string whitespace = " \\t\\r\\n";
    const std::size_t first = value.find_first_not_of( whitespace );
    if ( first == std::string::npos ) return "";
    const std::size_t last = value.find_last_not_of( whitespace );
    return value.substr( first, last - first + 1 );
}

std::string RemoveComment( const std::string & line )
{
    std::size_t end = line.size();
    const std::size_t hash = line.find( '#' );
    const std::size_t slash = line.find( "//" );
    if ( hash != std::string::npos ) end = std::min( end, hash );
    if ( slash != std::string::npos ) end = std::min( end, slash );
    return Trim( line.substr( 0, end ) );
}

bool IsParameterFileName( const std::string & name )
{
    const std::filesystem::path path( name );
    return ! name.empty()
        && path.filename() == path
        && name.find_first_of( "/\\\\:" ) == std::string::npos
        && ! path.is_absolute()
        && name != "." && name != ".."
        && path.extension() == ".txt"
        && name != "control.txt";
}

bool ParseManifestEntry( const std::string & source, std::string & name )
{
    const std::string line = RemoveComment( source );
    if ( line.empty() ) return true;

    std::string token = line;
    if ( line.front() == '"' || line.front() == '\'' )
    {
        const char quote = line.front();
        const std::size_t closing = line.find( quote, 1 );
        if ( closing == std::string::npos ) return false;
        name = line.substr( 1, closing - 1 );
        const std::string suffix = Trim( line.substr( closing + 1 ) );
        if ( ! suffix.empty() && suffix != ";" ) return false;
    }
    else
    {
        if ( token.back() == ';' )
            token = Trim( token.substr( 0, token.size() - 1 ) );
        if ( token.find_first_of( " \\t=(){}[]" ) != std::string::npos )
            return false;
        name = token;
    }
    return IsParameterFileName( name );
}

bool HasParameterDeclaration( const std::filesystem::path & path )
{
    std::ifstream input( path );
    if ( ! input ) return true;
    std::string line;
    while ( std::getline( input, line ) )
    {
        if ( RemoveComment( line ).find( '=' ) != std::string::npos )
            return true;
    }
    return false;
}

bool ReadFilenameManifest(
    const std::filesystem::path & controlPath,
    std::vector< std::string > & names )
{
    std::ifstream input( controlPath );
    if ( ! input ) return false;
    std::string line;
    bool firstLine = true;
    while ( std::getline( input, line ) )
    {
        if ( firstLine && line.size() >= 3
            && static_cast<unsigned char>( line[ 0 ] ) == 0xEF
            && static_cast<unsigned char>( line[ 1 ] ) == 0xBB
            && static_cast<unsigned char>( line[ 2 ] ) == 0xBF )
        {
            line.erase( 0, 3 );
        }
        firstLine = false;
        std::string name;
        if ( ! ParseManifestEntry( line, name ) ) return false;
        if ( ! name.empty() ) names.push_back( name );
    }
    return ! names.empty();
}

std::vector< std::filesystem::path > ListParameterFiles(
    const std::filesystem::path & scriptDirectory,
    bool includeControl )
{
    std::vector< std::filesystem::path > paths;
    std::error_code error;
    if ( ! std::filesystem::is_directory( scriptDirectory, error ) || error )
        return paths;

    for ( std::filesystem::directory_iterator it( scriptDirectory, error ), end;
          ! error && it != end; it.increment( error ) )
    {
        std::error_code statusError;
        if ( ! it->is_regular_file( statusError ) || statusError ) continue;
        const std::filesystem::path path = it->path();
        if ( path.extension() != ".txt" ) continue;
        if ( ! includeControl && path.filename() == "control.txt" ) continue;
        paths.push_back( path );
    }

    std::sort( paths.begin(), paths.end(),
        []( const std::filesystem::path & left,
            const std::filesystem::path & right )
        {
            return left.filename().generic_string()
                < right.filename().generic_string();
        } );
    return paths;
}

} // namespace

ScriptFileDiscovery DiscoverCaseScriptFiles( const std::string & caseDir )
{
    ScriptFileDiscovery result;
    const std::filesystem::path scriptDirectory =
        std::filesystem::path( caseDir ) / "script";
    const std::filesystem::path controlPath = scriptDirectory / "control.txt";

    std::error_code error;
    const bool hasControl =
        std::filesystem::is_regular_file( controlPath, error ) && ! error;

    if ( ! hasControl )
    {
        result.mode = ScriptFileSourceMode::Directory;
        result.diagnostics.push_back(
            "script/control.txt is absent; reading parameter files discovered in script/." );
        for ( const auto & path : ListParameterFiles( scriptDirectory, true ) )
        {
            if ( HasParameterDeclaration( path ) )
                result.filePaths.push_back( path.lexically_normal().string() );
            else
                result.diagnostics.push_back(
                    "Auxiliary script file '" + path.filename().string()
                    + "' has no parameter declarations and is not loaded as a parameter file." );
        }
        return result;
    }

    std::vector< std::string > manifestNames;
    if ( ! ReadFilenameManifest( controlPath, manifestNames ) )
    {
        result.mode = ScriptFileSourceMode::ParameterFile;
        result.diagnostics.push_back(
            "script/control.txt is not a filename-only manifest; reading it as a parameter file." );
        for ( const auto & path : ListParameterFiles( scriptDirectory, true ) )
        {
            if ( path.filename() == "control.txt" || HasParameterDeclaration( path ) )
                result.filePaths.push_back( path.lexically_normal().string() );
            else
                result.diagnostics.push_back(
                    "Auxiliary script file '" + path.filename().string()
                    + "' has no parameter declarations and is not loaded as a parameter file." );
        }
        return result;
    }

    result.mode = ScriptFileSourceMode::Manifest;
    const auto actualPaths = ListParameterFiles( scriptDirectory, false );
    std::set< std::string > actualNames;
    std::set< std::string > allNames;
    for ( const auto & path : actualPaths )
    {
        const std::string name = path.filename().generic_string();
        allNames.insert( name );
        if ( HasParameterDeclaration( path ) )
            actualNames.insert( name );
        else
            result.diagnostics.push_back(
                "Auxiliary script file '" + name
                + "' has no parameter declarations and is not loaded as a parameter file." );
    }

    std::set< std::string > declaredNames;
    for ( const std::string & name : manifestNames )
    {
        const bool firstOccurrence = declaredNames.insert( name ).second;
        if ( ! firstOccurrence )
        {
            result.diagnostics.push_back(
                "script/control.txt lists '" + name + "' more than once." );
            continue;
        }
        if ( actualNames.find( name ) == actualNames.end() )
        {
            if ( allNames.find( name ) != allNames.end() )
                result.diagnostics.push_back(
                    "script/control.txt lists '" + name
                    + "', but it has no parameter declarations and is not loaded." );
            else
                result.diagnostics.push_back(
                    "script/control.txt lists missing parameter file '" + name + "'." );
            continue;
        }
        result.filePaths.push_back(
            ( scriptDirectory / name ).lexically_normal().string() );
    }

    for ( const auto & path : actualPaths )
    {
        const std::string name = path.filename().generic_string();
        if ( actualNames.find( name ) != actualNames.end()
            && declaredNames.find( name ) == declaredNames.end() )
        {
            result.diagnostics.push_back(
                "Parameter file '" + name
                + "' exists in script/ but is not listed in control.txt; it will still be read." );
            result.filePaths.push_back( path.lexically_normal().string() );
        }
    }
    return result;
}

void ReadScriptFileNameList( std::vector< std::string > & scriptFileNameList )
{
    const ScriptFileDiscovery discovery = DiscoverCaseScriptFiles( Prj::prjBaseDir );
    scriptFileNameList.insert( scriptFileNameList.end(),
        discovery.filePaths.begin(), discovery.filePaths.end() );
}

void ReadScriptFileNameList(
    const std::string & caseDir,
    std::vector< std::string > & scriptFileNameList )
{
    const ScriptFileDiscovery discovery = DiscoverCaseScriptFiles( caseDir );
    scriptFileNameList.insert( scriptFileNameList.end(),
        discovery.filePaths.begin(), discovery.filePaths.end() );
}

} // namespace ONEFLOW

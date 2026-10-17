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

#include "Prj.h"
#include "Fatal.h"
#include "OStream.h"
#include "FileUtils.h"
#include <iostream>
#include <filesystem>
#include <stdexcept>

BeginNameSpace( ONEFLOW )

namespace
{

std::string NormalizePathString( std::filesystem::path path )
{
    path = path.lexically_normal();
    path.make_preferred();
    return path.string();
}

std::string FormatDirectoryPath( const std::string & directory )
{
    if ( directory.empty() )
    {
        return "";
    }

    std::filesystem::path path( directory );
    path = path.lexically_normal();
    path.make_preferred();

    std::string value = path.string();
    std::filesystem::path rootPath = path.root_path();
    rootPath.make_preferred();
    const std::string root = rootPath.string();

    // Remove redundant trailing separators, but preserve filesystem roots.
    while ( value.size() > root.size()
        && ! value.empty()
        && ( value.back() == '/' || value.back() == '\\' ) )
    {
        value.pop_back();
    }

    return value;
}

std::string ResolvePathUnderRoot(
    const std::string & rootName,
    const std::string & fileName )
{
    const std::filesystem::path path( fileName );

    // Only fully absolute paths bypass the configured root.
    if ( path.is_absolute() )
    {
        return NormalizePathString( path );
    }

    return NormalizePathString( std::filesystem::path( rootName ) / path );
}

}


bool Prj::hx_debug = false;
bool Prj::run_from_ide = false;
std::string Prj::system_root = "";
std::string Prj::execute_dir = "";
std::string Prj::current_dir = "";
std::string Prj::prjBaseDir = "";
std::string Prj::inputBaseDir = "";
std::string Prj::outputBaseDir = "";
std::string Prj::projectBaseDir = "";
std::string Prj::meshInputBaseDir = "";
std::string Prj::flowInputBaseDir = "";
std::string Prj::meshOutputBaseDir = "";
std::string Prj::flowOutputBaseDir = "";

void Prj::ProcessCmdLineArgs( std::vector<std::string> & args )
{
    CmdLineOptions opt = Prj::ParseCmdLineArgs( args );

    Prj::hx_debug = opt.debug;
    Prj::run_from_ide = opt.debug;

    // Process-level path initialization is performed once by Simulation.
    Prj::SetPrjBaseDir( opt.caseDir );
}

bool Prj::IsSystemRoot( const std::filesystem::path & path )
{
    std::error_code ec;

    if ( ! std::filesystem::is_directory( path, ec ) || ec )
    {
        return false;
    }

    const std::filesystem::path actionFile =
        path / "action" / "actionFileList.txt";

    ec.clear();

    return std::filesystem::is_regular_file( actionFile, ec ) && ! ec;
}

std::string Prj::FindSystemRoot()
{
    std::error_code ec;

    std::filesystem::path path( Prj::execute_dir );

    if ( path.empty() )
    {
        return "";
    }

    path = std::filesystem::absolute( path, ec );

    if ( ec )
    {
        return "";
    }

    while ( true )
    {
        std::filesystem::path candidate = path / "system";

        if ( Prj::IsSystemRoot( candidate ) )
        {
            std::string systemRoot = candidate.lexically_normal().string();

            if ( ! EndWithSlash( systemRoot ) )
            {
                systemRoot += "/";
            }

            return systemRoot;
        }

        std::filesystem::path parent = path.parent_path();

        if ( parent == path )
        {
            break;
        }

        path = parent;
    }

    return "";
}

void Prj::Init()
{
    Prj::execute_dir = HX_GetExeDirectory();
    Prj::current_dir = HX_GetCurrentDirectory();

    std::cout << " Prj::execute_dir = " << Prj::execute_dir << "\n";
    std::cout << " Prj::current_dir = " << Prj::current_dir << "\n";

    Prj::system_root = Prj::FindSystemRoot();

    if ( Prj::system_root.empty() )
    {
        Fatal( "Could not locate OneFLOW system directory from executable directory: "
            + Prj::execute_dir );
    }

    std::cout << " Prj::system_root = " << Prj::system_root << "\n";
}

std::string Prj::ResolveCaseDir( const std::string & caseDir )
{
    std::filesystem::path projectPath( caseDir );

    if ( projectPath.empty() )
    {
        Fatal( "Project path cannot be empty." );
    }

    if ( projectPath.is_relative() )
    {
        if ( Prj::current_dir.empty() )
        {
            Fatal( "Current directory is not initialized; cannot resolve "
                "relative project path: " + caseDir );
        }

        projectPath =
            std::filesystem::path( Prj::current_dir ) / projectPath;
    }

    projectPath = projectPath.lexically_normal();

    std::string resolvedDir = projectPath.string();

    // Keep the trailing slash because existing IO code relies on it.
    if ( ! EndWithSlash( resolvedDir ) )
    {
        resolvedDir += "/";
    }

    return resolvedDir;
}

void Prj::SetPrjBaseDir( const std::string & prjName )
{
    Prj::prjBaseDir = Prj::ResolveCaseDir( prjName );
    Prj::inputBaseDir = Prj::prjBaseDir;
    Prj::outputBaseDir = Prj::prjBaseDir;
    Prj::projectBaseDir.clear();
    Prj::meshInputBaseDir.clear();
    Prj::flowInputBaseDir.clear();
    Prj::meshOutputBaseDir.clear();
    Prj::flowOutputBaseDir.clear();

    std::cout << " Prj::prjBaseDir = "
        << Prj::prjBaseDir << "\n";
}

void Prj::ClearPrjBaseDir()
{
    // A completed case must not leave its directory bound to the process.
    Prj::prjBaseDir.clear();
    Prj::inputBaseDir.clear();
    Prj::outputBaseDir.clear();
    Prj::projectBaseDir.clear();
    Prj::meshInputBaseDir.clear();
    Prj::flowInputBaseDir.clear();
    Prj::meshOutputBaseDir.clear();
    Prj::flowOutputBaseDir.clear();
}

namespace
{
std::string ResolveBaseDirectory( const std::string & rootSpec )
{
    if ( rootSpec.empty() || rootSpec == "default" || rootSpec == "case" )
    {
        return Prj::prjBaseDir;
    }

    if ( rootSpec == "project" )
    {
        if ( Prj::projectBaseDir.empty() )
        {
            throw std::logic_error( "Project input/output root has not been configured" );
        }
        return Prj::projectBaseDir;
    }

    const std::filesystem::path rootPath( rootSpec );
    if ( rootPath.is_absolute() )
    {
        return rootPath.lexically_normal().string();
    }

    return ( std::filesystem::path( Prj::prjBaseDir ) / rootPath )
        .lexically_normal().string();
}
}

void Prj::SetProjectBaseDir( const std::string & rootDir )
{
    if ( rootDir.empty() )
    {
        throw std::invalid_argument( "Project root cannot be empty" );
    }

    std::filesystem::path path( rootDir );
    if ( path.is_relative() )
    {
        if ( Prj::current_dir.empty() )
        {
            throw std::logic_error( "Current directory is not initialized; cannot resolve project root" );
        }
        path = std::filesystem::path( Prj::current_dir ) / path;
    }

    Prj::projectBaseDir = path.lexically_normal().string();
}

void Prj::SetIOBaseDirs( const std::string & inputRoot, const std::string & outputRoot )
{
    Prj::inputBaseDir = ResolveBaseDirectory( inputRoot );
    Prj::outputBaseDir = ResolveBaseDirectory( outputRoot );
    Prj::meshInputBaseDir = Prj::inputBaseDir;
    Prj::flowInputBaseDir = Prj::inputBaseDir;
    Prj::meshOutputBaseDir = Prj::outputBaseDir;
    Prj::flowOutputBaseDir = Prj::outputBaseDir;
}

void Prj::SetDomainIOBaseDirs(
    const std::string & meshInputRoot,
    const std::string & flowInputRoot,
    const std::string & meshOutputRoot,
    const std::string & flowOutputRoot )
{
    Prj::meshInputBaseDir = meshInputRoot.empty() || meshInputRoot == "default"
        ? Prj::GetInputBaseDir() : ResolveBaseDirectory( meshInputRoot );
    Prj::flowInputBaseDir = flowInputRoot.empty() || flowInputRoot == "default"
        ? Prj::GetInputBaseDir() : ResolveBaseDirectory( flowInputRoot );
    Prj::meshOutputBaseDir = meshOutputRoot.empty() || meshOutputRoot == "default"
        ? Prj::GetOutputBaseDir() : ResolveBaseDirectory( meshOutputRoot );
    Prj::flowOutputBaseDir = flowOutputRoot.empty() || flowOutputRoot == "default"
        ? Prj::GetOutputBaseDir() : ResolveBaseDirectory( flowOutputRoot );
}

std::string Prj::GetInputBaseDir()
{
    return Prj::inputBaseDir.empty() ? Prj::prjBaseDir : Prj::inputBaseDir;
}

std::string Prj::GetOutputBaseDir()
{
    return Prj::outputBaseDir.empty() ? Prj::prjBaseDir : Prj::outputBaseDir;
}

std::string Prj::GetMeshInputBaseDir()
{
    return Prj::meshInputBaseDir.empty() ? Prj::GetInputBaseDir() : Prj::meshInputBaseDir;
}

std::string Prj::GetFlowInputBaseDir()
{
    return Prj::flowInputBaseDir.empty() ? Prj::GetInputBaseDir() : Prj::flowInputBaseDir;
}

std::string Prj::GetMeshOutputBaseDir()
{
    return Prj::meshOutputBaseDir.empty() ? Prj::GetOutputBaseDir() : Prj::meshOutputBaseDir;
}

std::string Prj::GetFlowOutputBaseDir()
{
    return Prj::flowOutputBaseDir.empty() ? Prj::GetOutputBaseDir() : Prj::flowOutputBaseDir;
}

std::string Prj::GetInputFileName( const std::string & fileName )
{
    return ResolvePathUnderRoot( Prj::GetInputBaseDir(), fileName );
}

std::string Prj::GetOutputFileName( const std::string & fileName )
{
    return ResolvePathUnderRoot( Prj::GetOutputBaseDir(), fileName );
}

std::string Prj::GetMeshInputFileName( const std::string & fileName )
{
    return ResolvePathUnderRoot( Prj::GetMeshInputBaseDir(), fileName );
}

std::string Prj::GetFlowInputFileName( const std::string & fileName )
{
    return ResolvePathUnderRoot( Prj::GetFlowInputBaseDir(), fileName );
}

std::string Prj::GetMeshOutputFileName( const std::string & fileName )
{
    return ResolvePathUnderRoot( Prj::GetMeshOutputBaseDir(), fileName );
}

std::string Prj::GetFlowOutputFileName( const std::string & fileName )
{
    return ResolvePathUnderRoot( Prj::GetFlowOutputBaseDir(), fileName );
}

void Prj::WriteIOPathInfo( const std::string & fileName )
{
    std::fstream file;
    Prj::OpenPrjFile(
        file,
        fileName,
        std::ios_base::out | std::ios_base::trunc );

    file << "case_root=" << FormatDirectoryPath( Prj::prjBaseDir ) << "\n";
    file << "project_root=" << FormatDirectoryPath( Prj::projectBaseDir ) << "\n";
    file << "input_root=" << FormatDirectoryPath( Prj::GetInputBaseDir() ) << "\n";
    file << "output_root=" << FormatDirectoryPath( Prj::GetOutputBaseDir() ) << "\n";
    file << "mesh_input_root=" << FormatDirectoryPath( Prj::GetMeshInputBaseDir() ) << "\n";
    file << "flow_input_root=" << FormatDirectoryPath( Prj::GetFlowInputBaseDir() ) << "\n";
    file << "mesh_output_root=" << FormatDirectoryPath( Prj::GetMeshOutputBaseDir() ) << "\n";
    file << "flow_output_root=" << FormatDirectoryPath( Prj::GetFlowOutputBaseDir() ) << "\n";
    file << "path_rule=absolute file paths bypass configured roots\n";

    Prj::CloseOutputFile( file, "WriteIOPathInfo" );
}

void Prj::OpenCaseFile(
    std::fstream & file,
    const std::string & caseDir,
    const std::string & fileName,
    const std::ios_base::openmode & openMode )
{
    std::string caseFileName = Prj::GetCaseFileName( caseDir, fileName );
    Prj::OpenFileWithParentDirectory( file, caseFileName, openMode );
}

void Prj::OpenPrjFile(
    std::fstream & file,
    const std::string & fileName,
    const std::ios_base::openmode & openMode )
{
    const bool writesOutput = ( openMode & std::ios_base::out ) != 0;
    const std::string resolvedFileName = writesOutput
        ? Prj::GetOutputFileName( fileName )
        : Prj::GetInputFileName( fileName );
    Prj::OpenFileWithParentDirectory( file, resolvedFileName, openMode );
}

void Prj::OpenFileWithParentDirectory(
    std::fstream & file,
    const std::string & fileName,
    const std::ios_base::openmode & openMode )
{
    // Create parent directories only for write operations.
    if ( ( openMode & std::ios_base::out ) != 0 )
    {
        CreateDirIfNeeded( fileName );
    }

    Prj::OpenFile( file, fileName, openMode );
}

void Prj::OpenFile(
    std::fstream & file,
    const std::string & fileName,
    const std::ios_base::openmode & openMode )
{
    file.open( fileName.c_str(), openMode );

    if ( ! file )
    {
        Fatal( "could not open " + fileName );
    }
}

void Prj::CloseFile( std::fstream & file )
{
    file.close();
    file.clear();
}

void Prj::CloseOutputFile( std::fstream & file, const std::string & operation )
{
    file.flush();
    if ( ! file )
    {
        throw std::runtime_error( operation + ": failed to write output data" );
    }

    file.close();
    if ( ! file )
    {
        throw std::runtime_error( operation + ": failed to close output file" );
    }
}

void Prj::MakePrjDir( const std::string & dirName )
{
    std::string prjDirName = Prj::GetPrjFileName( dirName );

    if ( ! HX_CreateDirectory( prjDirName ) )
    {
        Fatal( "Could not create project directory: " + prjDirName );
    }
}

// Same pattern as GetPrjFileName, but rooted at the OneFLOW installation's
// system directory (Prj::system_root) instead of the current case directory.
// Centralizing this here removes the scattered "Prj::system_root + ..."
// string concatenation that used to live in individual business-logic files.
std::string Prj::GetSystemFileName( const std::string & fileName )
{
    return ResolvePathUnderRoot( Prj::system_root, fileName );
}

std::string Prj::GetDirName( const std::string & fileName )
{
    size_t pos = fileName.find_last_of( "\\/" );

    if ( std::string::npos == pos )
    {
        return "";
    }
    else
    {
        return fileName.substr( 0, pos );
    }
}

void Prj::CreateDirIfNeeded( const std::string & prjFileName )
{
    std::string dirName = Prj::GetDirName( prjFileName );

    // A file in the current directory has no parent directory to create.
    if ( dirName.empty() )
    {
        return;
    }

    if ( ! HX_CreateDirectory( dirName ) )
    {
        Fatal( "Could not create parent directory for file: " + prjFileName );
    }
}

std::string Prj::GetPrjFileName( const std::string & fileName )
{
    return ResolvePathUnderRoot( Prj::prjBaseDir, fileName );
}

std::string Prj::GetCaseFileName(
    const std::string & caseDir,
    const std::string & fileName )
{
    std::filesystem::path baseDir( caseDir );

    if ( baseDir.is_relative() )
    {
        // A relative external case is relative to the currently active case,
        // not to the process working directory.
        baseDir =
            std::filesystem::path( Prj::prjBaseDir ) / baseDir;
    }

    return ResolvePathUnderRoot( baseDir.string(), fileName );
}

EndNameSpace
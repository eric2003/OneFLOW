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

BeginNameSpace( ONEFLOW )

bool Prj::hx_debug = false;
bool Prj::run_from_ide = false;
std::string Prj::system_root = "";
std::string Prj::execute_dir = "";
std::string Prj::current_dir = "";
std::string Prj::prjBaseDir = "";

Prj::Prj()
{
    ;
}

Prj::~Prj()
{
    ;
}

void Prj::ProcessCmdLineArgs( std::vector<std::string> & args )
{
    CmdLineOptions opt = Prj::ParseCmdLineArgs( args );

    if ( opt.debug )
    {
        Prj::hx_debug = true;
        Prj::run_from_ide = true;
    }

    Prj::Init();
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

void Prj::SetPrjBaseDir( const std::string & prjName )
{
    std::filesystem::path projectPath( prjName );

    if ( projectPath.empty() )
    {
        Fatal( "Project path cannot be empty." );
    }

    if ( projectPath.is_relative() )
    {
        projectPath =
            std::filesystem::path( Prj::current_dir ) / projectPath;
    }

    projectPath = projectPath.lexically_normal();

    Prj::prjBaseDir = projectPath.string();

    // Keep the trailing slash because existing IO code relies on it.
    if ( ! EndWithSlash( Prj::prjBaseDir ) )
    {
        Prj::prjBaseDir += "/";
    }

    std::cout << " Prj::prjBaseDir = "
        << Prj::prjBaseDir << "\n";
}

void Prj::OpenPrjFile(
    std::fstream & file,
    const std::string & fileName,
    const std::ios_base::openmode & openMode )
{
    std::string prjFileName = Prj::GetPrjFileName( fileName );

    // Create parent directories only for write operations.
    if ( ( openMode & std::ios_base::out ) != 0 )
    {
        CreateDirIfNeeded( prjFileName );
    }

    Prj::OpenFile( file, prjFileName, openMode );
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

void Prj::MakePrjDir( const std::string & dirName )
{
    std::string prjDirName = Prj::GetPrjFileName( dirName );

    HX_CreateDirectory( prjDirName );
}

// Same pattern as GetPrjFileName, but rooted at the OneFLOW installation's
// system directory (Prj::system_root) instead of the current case directory.
// Centralizing this here removes the scattered "Prj::system_root + ..."
// string concatenation that used to live in individual business-logic files.
std::string Prj::GetSystemFileName( const std::string & fileName )
{
    std::string fileNameNew = RemoveFirstSlash( fileName );

    return Prj::system_root + fileNameNew;
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

    if ( ! HX_IsDirectory( dirName ) )
    {
        HX_CreateDirectory( dirName );
    }
}

std::string Prj::GetPrjFileName( const std::string & fileName )
{
    std::string fileNameNew = RemoveFirstSlash( fileName );

    return Prj::prjBaseDir + fileNameNew;
}

EndNameSpace
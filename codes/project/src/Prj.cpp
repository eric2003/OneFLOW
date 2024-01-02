/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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
#include "Stop.h"
#include "OStream.h"
#include "FileUtil.h"
#include <iostream>

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

void Prj::ProcessCmdLineArgs( std::vector<std::string> &args )
{
    std::string choise = args[ 1 ];
    std::string prjName = args[ 2 ];
    if ( choise == "d" )
    {
        Prj::hx_debug = true;
        Prj::run_from_ide = true;
    }
    Prj::Init();
    Prj::SetPrjBaseDir( prjName );
}

void Prj::Init()
{
    Prj::execute_dir = HX_GetExePath();
    Prj::current_dir = HX_GetCurrentDir();

    std::cout << " Prj::execute_dir = " << Prj::execute_dir << "\n";
    std::cout << " Prj::current_dir = " << Prj::current_dir << "\n";

    std::string local_root = "/system/";
    if ( Prj::run_from_ide )
    {
        std::string current_dir_now = RemoveEndSlash( Prj::current_dir );
        Prj::system_root = current_dir_now + local_root;
    }
    else
    {
        std::string execute_dir = RemoveEndSlash( Prj::execute_dir );
        Prj::system_root = Prj::execute_dir + local_root;
    }
    std::cout << " Prj::system_root = " << Prj::system_root << "\n";
}

void Prj::SetPrjBaseDir( const std::string & prjName )
{
    std::string current_dir_now = RemoveEndSlash( Prj::current_dir );
    std::string prj_name_now = RemoveFirstSlash( prjName );
    ONEFLOW::StrIO << current_dir_now << "/" << prj_name_now;
    if ( ! EndWithSlash( prj_name_now ) )
    {
        ONEFLOW::StrIO << "/";
    }
    Prj::prjBaseDir = ONEFLOW::StrIO.str();
    std::cout << " Prj::prjBaseDir = " << Prj::prjBaseDir << "\n";
}

void Prj::OpenPrjFile( std::fstream & file, const std::string & fileName, const std::ios_base::openmode & openMode )
{
    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << Prj::prjBaseDir << fileName;

    std::string prjFileName = ONEFLOW::StrIO.str();

    CreateDirIfNeeded( prjFileName );

    Prj::OpenFile( file, prjFileName, openMode );
}

void Prj::OpenFile( std::fstream & file, const std::string & fileName, const std::ios_base::openmode & openMode )
{
    file.open( fileName.c_str(), openMode );
    if ( ! file )
    {
        std::cout << "could not open " << fileName << std::endl;
        Stop( "" );
    }
}

void Prj::CloseFile( std::fstream & file )
{
    file.close();
    file.clear();
}

void Prj::MakePrjDir( const std::string & dirName )
{
    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << Prj::prjBaseDir << dirName;

    std::string prjDirName = ONEFLOW::StrIO.str();
    //std::cout << " prjDirName = " << prjDirName << "\n";

    MakeDir( prjDirName );
}

std::string Prj::GetPrjDirName( const std::string & fileName )
{
    size_t pos = fileName.find_last_of("\\/");
    if ( std::string::npos == pos )
    {
        return "";
    }
    else
    {
        return fileName.substr(0, pos);
    }
}


void Prj::CreateDirIfNeeded( std::string & prjFileName )
{
    std::string prj_dir = Prj::GetPrjDirName( prjFileName );

    if ( ! DirExist( prj_dir ) )
    {
        MakeDir( prj_dir );
    }
}

std::string Prj::GetPrjFileName( const std::string & fileName )
{
    ONEFLOW::StrIO.ClearAll();

    std::string fileNameNew = RemoveFirstSlash( fileName );

    ONEFLOW::StrIO << Prj::prjBaseDir << fileNameNew;

    std::string prjFileName = ONEFLOW::StrIO.str();

    return prjFileName;
}

EndNameSpace

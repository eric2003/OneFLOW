/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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
#include "SimuCtrl.h"
#include "FileUtil.h"
#include <iostream>

BeginNameSpace( ONEFLOW )

std::string PrjStatus::prjBaseDir = "";

PrjStatus::PrjStatus()
{
    ;
}

PrjStatus::~PrjStatus()
{
    ;
}

void PrjStatus::SetPrjBaseDir( const std::string & prjName )
{
    std::string current_dir_now = RemoveEndSlash( SimuCtrl::current_dir );
    std::string prj_name_now = RemoveFirstSlash( prjName );
    ONEFLOW::StrIO << current_dir_now << "/" << prj_name_now;
    if ( ! EndWithSlash( prj_name_now ) )
    {
        ONEFLOW::StrIO << "/";
    }
    PrjStatus::prjBaseDir = ONEFLOW::StrIO.str();
    std::cout << " PrjStatus::prjBaseDir = " << PrjStatus::prjBaseDir << "\n";
}

void MakePrjDir( const std::string & dirName )
{
    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << PrjStatus::prjBaseDir << dirName;

    std::string prjDirName = ONEFLOW::StrIO.str();
    //cout << " prjDirName = " << prjDirName << "\n";

    MakeDir( prjDirName );
}

std::string GetPrjDirName( const std::string & fileName )
{
    size_t pos = fileName.find_last_of("\\/");
    if ( string::npos == pos )
    {
        return "";
    }
    else
    {
        return fileName.substr(0, pos);
    }
}


void CreateDirIfNeeded( std::string & prjFileName )
{
    std::string prj_dir = ONEFLOW::GetPrjDirName( prjFileName );

    if ( ! DirExist( prj_dir ) )
    {
        MakeDir( prj_dir );
    }
}

void OpenPrjFile( std::fstream & file, const std::string & fileName, const std::ios_base::openmode & openMode )
{
    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << PrjStatus::prjBaseDir << fileName;

    std::string prjFileName = ONEFLOW::StrIO.str();

    CreateDirIfNeeded( prjFileName );

    ONEFLOW::OpenFile( file, prjFileName, openMode );
}

std::string GetPrjFileName( const std::string & fileName )
{
    ONEFLOW::StrIO.ClearAll();

    std::string fileNameNew = RemoveFirstSlash( fileName );

    ONEFLOW::StrIO << PrjStatus::prjBaseDir << fileNameNew;

    std::string prjFileName = ONEFLOW::StrIO.str();

    return prjFileName;
}

void OpenFile( std::fstream & file, const std::string & fileName, const std::ios_base::openmode & openMode )
{
    file.open( fileName.c_str(), openMode );
    if ( ! file )
    {
        std::cout << "could not open " << fileName << std::endl;
        Stop( "" );
    }
}

void CloseFile( std::fstream & file )
{
    file.close();
    file.clear();
}

EndNameSpace

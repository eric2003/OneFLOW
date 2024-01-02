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
#pragma once
#include "Configure.h"
#include <fstream>
#include <string>
#include <vector>

BeginNameSpace( ONEFLOW )

class Prj
{
public:
    Prj();
    ~Prj();
public:
    static bool hx_debug;
    static bool run_from_ide;
    static std::string system_root;
    static std::string current_dir;
    static std::string execute_dir;
    static std::string prjBaseDir;
public:
    static void Init();
    static void SetPrjBaseDir( const std::string & prjName );
    static void ProcessCmdLineArgs( std::vector<std::string> &args );
public:
    static void OpenPrjFile( std::fstream & file, const std::string & fileName, const std::ios_base::openmode & openMode );
    static void OpenFile( std::fstream & file, const std::string & fileName, const std::ios_base::openmode & openMode );
    static void CloseFile( std::fstream & file );
    static void CreateDirIfNeeded( std::string & prjFileName );
    static std::string GetPrjFileName( const std::string & fileName );
    static std::string GetPrjDirName( const std::string & fileName );
    static void MakePrjDir( const std::string & dirName );
};



EndNameSpace

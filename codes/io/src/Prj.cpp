/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
#include "OStream.h"
#include "FileUtil.h"
#include "LogFile.h"
#include "SimuCtrl.h"


#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

string PrjStatus::prjBaseDir = "";

PrjStatus::PrjStatus()
{
    ;
}

PrjStatus::~PrjStatus()
{
    ;
}

void PrjStatus::SetPrjBaseDir( const string & prjName )
{
    string current_dir_now = RemoveEndSlash( SimuCtrl::current_dir );
    string prj_name_now = RemoveFirstSlash( prjName );
    ONEFLOW::StrIO << current_dir_now << "/" << prj_name_now;
    if ( ! EndWithSlash( prj_name_now ) )
    {
        ONEFLOW::StrIO << "/";
    }
    PrjStatus::prjBaseDir = ONEFLOW::StrIO.str();
    cout << " PrjStatus::prjBaseDir =  " << PrjStatus::prjBaseDir << "\n";
}

void MakePrjDir( const string & dirName )
{
    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << PrjStatus::prjBaseDir << dirName;

    string prjDirName = ONEFLOW::StrIO.str();
    //cout << " prjDirName = " << prjDirName << "\n";

    MakeDir( prjDirName );
}

void OpenPrjFile( fstream & file, const string & fileName, const ios_base::openmode & openMode )
{
    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << PrjStatus::prjBaseDir << fileName;

    string prjFileName = ONEFLOW::StrIO.str();

    ONEFLOW::OpenFile( file, prjFileName, openMode );
}

string GetPrjFileName( const string & fileName )
{
    ONEFLOW::StrIO.ClearAll();

    string fileNameNew = RemoveFirstSlash( fileName );

    ONEFLOW::StrIO << PrjStatus::prjBaseDir << fileNameNew;

    string prjFileName = ONEFLOW::StrIO.str();

    return prjFileName;
}


EndNameSpace

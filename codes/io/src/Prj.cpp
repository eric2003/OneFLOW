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
    string baseDir = "./";
    ONEFLOW::StrIO << baseDir << prjName;
    if ( ! EndWithForwardSlash( prjName ) )
    {
        ONEFLOW::StrIO << "/";
    }
    PrjStatus::prjBaseDir = ONEFLOW::StrIO.str();
    cout << " PrjStatus::prjBaseDir =  " << PrjStatus::prjBaseDir << "\n";
}

bool EndWithForwardSlash( const string & fileName )
{
    size_t pos = fileName.find_last_of("/");
    size_t ss = fileName.size();
    if ( ss == 0 )
    {
        return false;
    }
    else
    {
        bool flag = fileName.substr( ss - 1, 1 ) == "/";
        return flag;
    }
}

bool StartWithForwardSlash( const string & fileName )
{
    size_t pos = fileName.find_first_of("/");
    if ( fileName.size() == 0 )
    {
        return false;
    }

    if ( fileName.substr( 0,1 ) == "/" )
    {
        return true;
    }
    return false;
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

string RemoveFirstSlash( const string & fileName )
{
    if ( StartWithForwardSlash( fileName ) )
    {
        int len = fileName.size();
        return fileName.substr( 1, len - 1 );
    }
    return fileName;
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

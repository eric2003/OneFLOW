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
#pragma once
#include "Configure.h"
#include <fstream>
#include <string>
#include <iomanip>
#include <sstream>
#include <set>
#include <vector>
#include <algorithm>
using namespace std;

BeginNameSpace( ONEFLOW )

bool DirExist( const string & dirName );
void MakeDir( const string & dirName );

string HX_GetExePath();
string HX_GetCurrentDir();

bool EndWithSlash( const string & fileName );
bool EndWithBackwardSlash( const string & fileName );
bool EndWithForwardSlash( const string & fileName );
bool StartWithForwardSlash( const string & fileName );
string RemoveFirstSlash( const string & fileName );
string RemoveEndSlash( const string & fileName );

void OpenFile( fstream & file, const string & fileName, const ios_base::openmode & openMode );
void CloseFile( fstream & file );

void   GetFileNameExtension( const string & fullName, string & mainName, string & extensionName, const string & fileNameSeparator );

void   ModifyFileMainName     ( string & fileName, const string & newMainName );
void   ModifyFileExtensionName( string & fileName, const string & newExtensionName );

template < typename T >
string AddSymbolToFileName( const string & fileName, const T & symbol )
{
    string mainName, extensionName;
    ONEFLOW::GetFileNameExtension( fileName, mainName, extensionName, "." );

    ostringstream oss;
    oss << mainName << symbol << "." << extensionName;
    return oss.str();
}

template < typename T1, typename T2 >
string AddSymbolToFileName( const string & fileName, const T1 & v1, const T2 & v2 )
{
    string mainName, extensionName;
    ONEFLOW::GetFileNameExtension( fileName, mainName, extensionName, "." );

    ostringstream oss;
    oss << mainName << v1 << v2 << "." << extensionName;
    return oss.str();
}

template < typename T1, typename T2, typename T3 >
string AddSymbolToFileName( const string & fileName, const T1 & v1, const T2 & v2, const T3 & v3 )
{
    string mainName, extensionName;
    ONEFLOW::GetFileNameExtension( fileName, mainName, extensionName, "." );

    ostringstream oss;
    oss << mainName << v1 << v2 << v3 << "." << extensionName;
    return oss.str();
}


EndNameSpace

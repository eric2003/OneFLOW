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
#include <iomanip>
#include <sstream>
#include <set>
#include <vector>
#include <algorithm>

BeginNameSpace( ONEFLOW )

bool DirExist( const std::string & dirName );
void MakeDir( const std::string & dirName );

std::string HX_GetExePath();
std::string HX_GetCurrentDir();

bool EndWithSlash( const std::string & fileName );
bool EndWithBackwardSlash( const std::string & fileName );
bool EndWithForwardSlash( const std::string & fileName );
bool StartWithForwardSlash( const std::string & fileName );
std::string RemoveFirstSlash( const std::string & fileName );
std::string RemoveEndSlash( const std::string & fileName );

void   GetFileNameExtension( const std::string & fullName, std::string & mainName, std::string & extensionName, const std::string & fileNameSeparator );

void   ModifyFileMainName     ( std::string & fileName, const std::string & newMainName );
void   ModifyFileExtensionName( std::string & fileName, const std::string & newExtensionName );

template < typename T >
std::string AddSymbolToFileName( const std::string & fileName, const T & symbol )
{
    std::string mainName, extensionName;
    ONEFLOW::GetFileNameExtension( fileName, mainName, extensionName, "." );

    std::ostringstream oss;
    oss << mainName << symbol << "." << extensionName;
    return oss.str();
}

template < typename T1, typename T2 >
std::string AddSymbolToFileName( const std::string & fileName, const T1 & v1, const T2 & v2 )
{
    std::string mainName, extensionName;
    ONEFLOW::GetFileNameExtension( fileName, mainName, extensionName, "." );

    std::ostringstream oss;
    oss << mainName << v1 << v2 << "." << extensionName;
    return oss.str();
}

template < typename T1, typename T2, typename T3 >
std::string AddSymbolToFileName( const std::string & fileName, const T1 & v1, const T2 & v2, const T3 & v3 )
{
    std::string mainName, extensionName;
    ONEFLOW::GetFileNameExtension( fileName, mainName, extensionName, "." );

    std::ostringstream oss;
    oss << mainName << v1 << v2 << v3 << "." << extensionName;
    return oss.str();
}


EndNameSpace

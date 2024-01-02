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

#include "FileUtil.h"
#include "Stop.h"

#ifdef _WINDOWS
#include <windows.h>
#include <direct.h>
#include <io.h>
#else
    #ifdef WIN_GNU
        #include <windows.h>
        #include <direct.h>
        #include <io.h>
    #else
        #include <sys/stat.h>
        #include <unistd.h>
    #endif
#endif

#include <iostream>

BeginNameSpace( ONEFLOW )

bool DirExist( const std::string & dirName )
{
#ifdef _WINDOWS
    bool flag = ( _access( dirName.c_str(), 0 ) == 0 );
    return flag;
#else
    bool flag = ( access( dirName.c_str(), 0 ) == 0 );
    return flag;
#endif
}

void MakeDir( const std::string & dirName )
{
    int flag;
#ifdef _WINDOWS
    flag = _mkdir( dirName.c_str() );
#else
    #ifdef WIN_GNU
        flag = _mkdir( dirName.c_str() );
    #else
        flag = mkdir( dirName.c_str(), S_IRWXU );
    #endif
#endif

    if ( flag == 0 )
    {
        std::cout << dirName << " directory has been created successfully !\n";
    }
}

std::string HX_GetExePath()
{
    char buffer[ FILENAME_MAX ] = { 0 };
#ifdef _WIN32
    GetModuleFileName( NULL, buffer, FILENAME_MAX );
#else
    std::size_t count = readlink( "/proc/self/exe", buffer, FILENAME_MAX );
#endif
    std::string::size_type pos = std::string( buffer ).find_last_of( "\\/" );
    return std::string( buffer ).substr( 0, pos);
}


std::string HX_GetCurrentDir()
{
#ifdef _WINDOWS
    char * cwd = _getcwd( 0, 0 );
#else
    char * cwd = getcwd( 0, 0 ); 
#endif
    std::string working_dir( cwd ) ;
    std::free( cwd ) ;
    return working_dir ;
}

bool EndWithSlash( const std::string & fileName )
{
    if ( EndWithForwardSlash( fileName ) ||
        EndWithBackwardSlash( fileName ) )
    {
        return true;
    }
    return false;
}

bool EndWithBackwardSlash( const std::string & fileName )
{
    size_t pos = fileName.find_last_of("\\");
    size_t ss = fileName.size();
    if ( ss == 0 )
    {
        return false;
    }
    else
    {
        bool flag = fileName.substr( ss - 1, 1 ) == "\\";
        return flag;
    }
}

bool EndWithForwardSlash( const std::string & fileName )
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

bool StartWithForwardSlash( const std::string & fileName )
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

std::string RemoveFirstSlash( const std::string & fileName )
{
    if ( StartWithForwardSlash( fileName ) )
    {
        int len = fileName.size();
        return fileName.substr( 1, len - 1 );
    }
    return fileName;
}

std::string RemoveEndSlash( const std::string & fileName )
{
    if ( EndWithSlash( fileName ) )
    {
        int len = fileName.size();
        return fileName.substr( 0, len - 1 );
    }
    return fileName;

}

void GetFileNameExtension( const std::string & fullName, std::string & mainName, std::string & extensionName, const std::string & fileNameSeparator )
{
    std::basic_string <char>::size_type index;

    index         = fullName.find_last_of( fileNameSeparator );
    mainName      = fullName.substr( 0, index );
    extensionName = fullName.substr( index+1, fullName.length() - index - 1 );
}

void ModifyFileMainName( std::string & fileName,  const std::string & newMainName )
{
    std::string mainName, extensionName;
    ONEFLOW::GetFileNameExtension( fileName, mainName, extensionName, "." );

    std::ostringstream oss;
    oss << newMainName << "." << extensionName;

    fileName = oss.str();
}

void ModifyFileExtensionName( std::string & fileName,  const std::string & newExtensionName )
{
    std::string mainName, extensionName;
    ONEFLOW::GetFileNameExtension( fileName, mainName, extensionName, "." );

    std::ostringstream oss;
    oss << mainName << "." << newExtensionName;

    fileName = oss.str();
}

EndNameSpace

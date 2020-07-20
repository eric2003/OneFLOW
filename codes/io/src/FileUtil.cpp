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

#include "FileUtil.h"
#include "Word.h"
#include "Stop.h"
#include "BasicParallel.h"
#include "LogFile.h"
#ifdef _WINDOWS
#include <windows.h>
#include <direct.h>
#include <io.h>
#else
#include <sys/stat.h>
#include <unistd.h>
#endif

#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

bool DirExist( const string & dirName )
{
#ifdef _WINDOWS
    bool flag = ( _access( dirName.c_str(), 0 ) == 0 );
    return flag;
#else
    bool flag = ( access( dirName.c_str(), 0 ) == 0 );
    return flag;
#endif
}

void MakeDir( const string & dirName )
{
    int flag;
#ifdef _WINDOWS
    flag = _mkdir( dirName.c_str() );
#else
    flag = mkdir( dirName.c_str(), S_IRWXU );
#endif

    if ( flag == 0 )
    {
        cout << dirName << " directory has been created successfully !\n";
    }
}

string HX_GetExePath()
{
    char buffer[ FILENAME_MAX ] = { 0 };
#ifdef _WIN32
    GetModuleFileName( NULL, buffer, FILENAME_MAX );
#else
    ssize_t count = readlink( "/proc/self/exe", buffer, FILENAME_MAX );
#endif
    string::size_type pos = string( buffer ).find_last_of( "\\/" );
    return string( buffer ).substr( 0, pos);
}

string HX_GetCurrentDir()
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

bool EndWithSlash( const string & fileName )
{
    if ( EndWithForwardSlash( fileName ) ||
        EndWithBackwardSlash( fileName ) )
    {
        return true;
    }
    return false;
}

bool EndWithBackwardSlash( const string & fileName )
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

string RemoveFirstSlash( const string & fileName )
{
    if ( StartWithForwardSlash( fileName ) )
    {
        int len = fileName.size();
        return fileName.substr( 1, len - 1 );
    }
    return fileName;
}

string RemoveEndSlash( const string & fileName )
{
    if ( EndWithSlash( fileName ) )
    {
        int len = fileName.size();
        return fileName.substr( 0, len - 1 );
    }
    return fileName;

}

void OpenFile( fstream & file, const string & fileName, const ios_base::openmode & openMode )
{
    file.open( fileName.c_str(), openMode );
    if ( ! file )
    {
        cout << "could not open " << fileName << endl;
        Stop( "" );
    }
}

void CloseFile( fstream & file )
{
    file.close();
    file.clear();
}


void GetFileNameExtension( const string & fullName, string & mainName, string & extensionName, const string & fileNameSeparator )
{
    basic_string <char>::size_type index;

    index         = fullName.find_last_of( fileNameSeparator );
    mainName      = fullName.substr( 0, index );
    extensionName = fullName.substr( index+1, fullName.length() - index - 1 );
}

void ModifyFileMainName( string & fileName,  const string & newMainName )
{
    string mainName, extensionName;
    ONEFLOW::GetFileNameExtension( fileName, mainName, extensionName, "." );

    ostringstream oss;
    oss << newMainName << "." << extensionName;

    fileName = oss.str();
}

void ModifyFileExtensionName( string & fileName,  const string & newExtensionName )
{
    string mainName, extensionName;
    ONEFLOW::GetFileNameExtension( fileName, mainName, extensionName, "." );

    ostringstream oss;
    oss << mainName << "." << newExtensionName;

    fileName = oss.str();
}

EndNameSpace

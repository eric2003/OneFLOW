/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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

#include "BasicIO.h"
#include "Word.h"
#include "Stop.h"
#include "BasicParallel.h"
#include "LogFile.h"
#ifdef _WINDOWS
#include <direct.h>
#include <io.h>
#else
#include <sys/stat.h>
#include <unistd.h>
#endif

#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

OStream StrIO;

void OStream::ClearAll()
{
    this->str("");
}

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


void Clear( ostringstream & oss )
{
    oss.clear( );
    oss.str( "" );
}


//void StopProgramFunction( const ostringstream & oss, const string & fileName, const int & fileLine, const string & dateName, const string & timeName )
//{
//    ONEFLOW::StopProgramFunction( oss.str(), fileName, fileLine, dateName, timeName );
//}
//
//void StopProgramFunction( const string & stopInformation, const string & fileName, const int & fileLine, const string & dateName, const string & timeName )
//{
//    ONEFLOW::HXFinalize();
//    cout << endl;
//    cout << "++++++++++++++++++Stop Information  +++++++++++++++++++++++++++++\n";
//    cout <<  stopInformation << endl;
//    cout << " The stop filename is : " << fileName << endl;
//    cout << " at line " << fileLine << endl;
//    cout << " Compiled On line " << dateName << " at " << timeName << endl;
//    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
//    exit( 0 );
//}


EndNameSpace

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

void TrimBlanks( string & source )
{
    string::size_type firstIndex, nextIndex;

    firstIndex = source.find_first_not_of( " " );
    nextIndex  = source.find_last_not_of( " " );

    source = source.substr( firstIndex, nextIndex - firstIndex + 1 );
}

void SkipLines( fstream & file, int numberOfLines )
{
    string line;
    for ( int iLine = 0; iLine < numberOfLines; ++ iLine )
    {
        ONEFLOW::ReadNextLine( file, line );
        if ( file.eof() ) return;
    }
}

void ReadNextLine( fstream & file, string & line )
{
    std::getline( file, line );
}

string FindNextWord( const string & source, string & word, const string & separator )
{
    string::size_type firstIndex, nextIndex, notFindIndex = - 1;
    string emptyString = "";
    firstIndex = source.find_first_not_of( separator, 0 );

    if ( firstIndex == notFindIndex )
    {
        word = emptyString;
        return emptyString;
    }
    nextIndex = source.find_first_of( separator, firstIndex );
    if ( nextIndex == notFindIndex )
    {
        word = source.substr( firstIndex );
        return emptyString;
    }
    else
    {
        word = source.substr( firstIndex, nextIndex - firstIndex );
        return source.substr( nextIndex );
    }
    return emptyString;
}

string FindNextWord( string & source, const string & separator )
{
    string::size_type firstIndex, nextIndex, notFindIndex = - 1;
    string emptyString = "";
    firstIndex = source.find_first_not_of( separator, 0 );

    if ( firstIndex == notFindIndex )
    {
        return emptyString;
    }

    nextIndex = source.find_first_of( separator, firstIndex );
    if ( nextIndex == notFindIndex )
    {
        string word = source.substr( firstIndex );
        source = emptyString;
        return word;
    }
    else
    {
        string word = source.substr( firstIndex, nextIndex - firstIndex );
        source = source.substr( nextIndex );
        return word;
    }
    return emptyString;
}

bool FindString( const string & source, const string & word )
{
    return source.find( word ) != string::npos;

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

bool ReadNextNonEmptyLine( fstream & file, string & line )
{
    bool isSpaceLine = true;

    do
    {
        ONEFLOW::ReadNextLine( file, line );

        for ( string::iterator iter = line.begin(); iter != line.end(); ++ iter )
        {
            if ( ! isspace( * iter ) )
            {
                isSpaceLine = false;
                break;
            }
        }
        if ( file.eof() )
        {
            if ( isSpaceLine )
            {
                return false;
            }
            else
            {
                return true;
            }
        }
    }     while ( isSpaceLine );
    
    return true;
}

void Clear( ostringstream & oss )
{
    oss.clear( );
    oss.str( "" );
}


void StopProgramFunction( const ostringstream & oss, const string & fileName, const int & fileLine, const string & dateName, const string & timeName )
{
    ONEFLOW::StopProgramFunction( oss.str(), fileName, fileLine, dateName, timeName );
}

void StopProgramFunction( const string & stopInformation, const string & fileName, const int & fileLine, const string & dateName, const string & timeName )
{
    ONEFLOW::HXFinalize();
    cout << endl;
    cout << "++++++++++++++++++Stop Information  +++++++++++++++++++++++++++++\n";
    cout <<  stopInformation << endl;
    cout << " The stop filename is : " << fileName << endl;
    cout << " at line " << fileLine << endl;
    cout << " Compiled On line " << dateName << " at " << timeName << endl;
    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
    exit( 0 );
}

void ToLowerCase( string & word )
{
    std::transform( word.begin(), word.end(), word.begin(), StringToLowerCase() );
}

void ToUpperCase( string & word )
{
    std::transform( word.begin(), word.end(), word.begin(), StringToUpperCase() );
}

//int GetIntegerDigitWidth( int integerData )
//{
//    return 1 + static_cast< int >( floor( log10( static_cast< double > ( ONEFLOW::MAX( 1, integerData ) ) ) ) );
//}

bool IsDigit( const string & word )
{
    for ( int i = 0; i < word.size(); ++ i )
    {
        if ( ! isdigit( word[ i ] ) )
        {
            return false;
        }
    }
    return true;
}

EndNameSpace

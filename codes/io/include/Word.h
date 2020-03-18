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
#include "HXDefine.h"
#include <fstream>
using namespace std;

BeginNameSpace( ONEFLOW )

class StringToLowerCase
{
public:
    char operator()( char val )
    {
        return tolower( val );
    }
};

class StringToUpperCase
{
public:
    char operator()( char val )
    {
        return toupper( val );
    }
};


//Character conversion to numeric value
template < typename T >
inline T StringToDigit( const string & str, ios_base & ( * f )( ios_base & ) = std::dec )
{
    T value;
    istringstream stream( str );
    stream >> f >> value;
    return value;
}

//Convert numeric to character
template < typename T >
inline string DigitToString( T t, ios_base & ( * f )( ios_base & ) = std::dec )
{
    ostringstream oss;
    oss << f << t;
    return oss.str();
}

//int GetIntegerDigitWidth( int integerData );

class Word
{
public:
    Word();
    ~Word();
public:
    static void ReadNextLine( fstream & file, string & line );
    static void SkipLines( fstream & file, int numberOfLines );
    static void TrimBlanks( string & source );
    static bool FindString( const string & source, const string & word );
    static bool IsEmptyLine  ( const string & line );
    static bool IsCommentLine( const string & line );
    static bool IsCommentLine( const string & line, StringField & comlist );
public:
    static string TMP_FindNextWord( const string & source, string & word, const string & separator );
    static string FindNextWord( string & source, const string & separator );
public:
    static void ToLowerCase( string & word );
    static void ToUpperCase( string & word );
public:
    static bool ReadNextNonEmptyLine( fstream & file, string & line );
    static bool IsDigit( const string & word );

};


EndNameSpace
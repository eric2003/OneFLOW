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
#include <sstream>
#include <string>
#include <vector>

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
inline T StringToDigit( const std::string & str, std::ios_base & ( * f )( std::ios_base & ) = std::dec )
{
    T value;
    std::istringstream stream( str );
    stream >> f >> value;
    return value;
}

//Convert numeric to character
template < typename T >
inline std::string DigitToString( T t, std::ios_base & ( * f )( std::ios_base & ) = std::dec )
{
    std::ostringstream oss;
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
    static void ReadNextLine( std::fstream & file, std::string & line );
    static void SkipLines( std::fstream & file, int numberOfLines );
    static void TrimBlanks( std::string & source );
    static bool FindString( const std::string & source, const std::string & word );
    static bool IsEmptyLine  ( const std::string & line );
    static bool IsCommentLine( const std::string & line );
    static bool IsCommentLine( const std::string & line, std::vector<std::string> & comlist );
public:
    static std::string TMP_FindNextWord( const std::string & source, std::string & word, const std::string & separator );
    static std::string FindNextWord( std::string & source, const std::string & separator );
public:
    static void ToLowerCase( std::string & word );
    static void ToUpperCase( std::string & word );
public:
    static bool ReadNextNonEmptyLine( std::fstream & file, std::string & line );
    static bool IsDigit( const std::string & word );

};


EndNameSpace

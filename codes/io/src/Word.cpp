/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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

#include "Word.h"
#include <algorithm>

BeginNameSpace( ONEFLOW )

Word::Word()
{
    ;
}

Word::~Word()
{
    ;
}

void Word::SkipLines( std::fstream & file, int numberOfLines )
{
    std::string line;
    for ( int iLine = 0; iLine < numberOfLines; ++ iLine )
    {
        Word::ReadNextLine( file, line );
        if ( file.eof() ) return;
    }
}

void Word::ReadNextLine( std::fstream & file, std::string & line )
{
    std::getline( file, line );
}

void Word::TrimBlanks( std::string & source )
{
    std::string::size_type firstIndex, nextIndex;

    firstIndex = source.find_first_not_of( " " );
    nextIndex  = source.find_last_not_of( " " );

    source = source.substr( firstIndex, nextIndex - firstIndex + 1 );
}

bool Word::FindString( const std::string & source, const std::string & word )
{
    return source.find( word ) != std::string::npos;
}

std::string Word::TMP_FindNextWord( const std::string & source, std::string & word, const std::string & separator )
{
    std::string::size_type firstIndex, nextIndex, notFindIndex = std::string::npos;
    std::string emptyString = "";
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

bool Word::IsEmptyLine( const std::string & line )
{
    if ( line == "" )
    {
        return true;
    }
    else
    {
        const std::string notReadableSeparator = " \r\n\t";
        std::string word;
        Word::TMP_FindNextWord( line, word, notReadableSeparator );
        return word == "";
    }
}

bool Word::IsCommentLine( const std::string & line )
{
    const std::string notReadableSeparator = " \r\n\t";
    std::string word;
    Word::TMP_FindNextWord( line, word, notReadableSeparator );
    return ( word.substr( 0, 1 ) == "#" ||
        word.substr( 0, 2 ) == "//" );
}

bool Word::IsCommentLine(const std::string& line, std::vector<std::string> &comlist )
{
    const std::string notReadableSeparator = " \r\n\t";
    std::string word;
    Word::TMP_FindNextWord(line, word, notReadableSeparator);
    for (int i = 0; i < comlist.size(); ++ i)
    {
        std::string & t = comlist[ i ];
        int n = t.length();
        if ( word.substr(0, n) == t ) return true;
    }
    return false;
}

std::string Word::FindNextWord( std::string & source, const std::string & separator )
{
    std::string::size_type firstIndex, nextIndex, notFindIndex = std::string::npos;
    std::string emptyString = "";
    firstIndex = source.find_first_not_of( separator, 0 );

    if ( firstIndex == notFindIndex )
    {
        return emptyString;
    }

    nextIndex = source.find_first_of( separator, firstIndex );
    if ( nextIndex == notFindIndex )
    {
        std::string word = source.substr( firstIndex );
        source = emptyString;
        return word;
    }
    else
    {
        std::string word = source.substr( firstIndex, nextIndex - firstIndex );
        source = source.substr( nextIndex );
        return word;
    }
    return emptyString;
}

void Word::ToLowerCase( std::string & word )
{
    std::transform( word.begin(), word.end(), word.begin(), StringToLowerCase() );
}

void Word::ToUpperCase( std::string & word )
{
    std::transform( word.begin(), word.end(), word.begin(), StringToUpperCase() );
}

bool Word::ReadNextNonEmptyLine( std::fstream & file, std::string & line )
{
    bool isSpaceLine = true;

    do
    {
        Word::ReadNextLine( file, line );

        for ( std::string::iterator iter = line.begin(); iter != line.end(); ++ iter )
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

bool Word::IsDigit( const std::string & word )
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

//int GetIntegerDigitWidth( int integerData )
//{
//    return 1 + static_cast< int >( floor( log10( static_cast< double > ( ONEFLOW::MAX( 1, integerData ) ) ) ) );
//}

EndNameSpace

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

#include "Word.h"
#include "BasicIO.h"

BeginNameSpace( ONEFLOW )

Word::Word()
{
    ;
}

Word::~Word()
{
    ;
}

void Word::SkipLines( fstream & file, int numberOfLines )
{
    string line;
    for ( int iLine = 0; iLine < numberOfLines; ++ iLine )
    {
        Word::ReadNextLine( file, line );
        if ( file.eof() ) return;
    }
}

void Word::ReadNextLine( fstream & file, string & line )
{
    std::getline( file, line );
}

void Word::TrimBlanks( string & source )
{
    string::size_type firstIndex, nextIndex;

    firstIndex = source.find_first_not_of( " " );
    nextIndex  = source.find_last_not_of( " " );

    source = source.substr( firstIndex, nextIndex - firstIndex + 1 );
}

bool Word::FindString( const string & source, const string & word )
{
    return source.find( word ) != string::npos;
}

string Word::TMP_FindNextWord( const string & source, string & word, const string & separator )
{
    string::size_type firstIndex, nextIndex, notFindIndex = string::npos;
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

bool Word::IsEmptyLine( const string & line )
{
    if ( line == "" )
    {
        return true;
    }
    else
    {
        const string notReadableSeparator = " \r\n\t";
        string word;
        Word::TMP_FindNextWord( line, word, notReadableSeparator );
        return word == "";
    }
}

bool Word::IsCommentLine( const string & line )
{
    const string notReadableSeparator = " \r\n\t";
    string word;
    Word::TMP_FindNextWord( line, word, notReadableSeparator );
    return ( word.substr( 0, 1 ) == "#" ||
        word.substr( 0, 2 ) == "//" );
}

bool Word::IsCommentLine(const string& line, StringField &comlist )
{
    const string notReadableSeparator = " \r\n\t";
    string word;
    Word::TMP_FindNextWord(line, word, notReadableSeparator);
    for (int i = 0; i < comlist.size(); ++ i)
    {
        string & t = comlist[ i ];
        int n = t.length();
        if ( word.substr(0, n) == t ) return true;
    }
    return false;
}

string Word::FindNextWord( string & source, const string & separator )
{
    string::size_type firstIndex, nextIndex, notFindIndex = string::npos;
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

void Word::ToLowerCase( string & word )
{
    std::transform( word.begin(), word.end(), word.begin(), StringToLowerCase() );
}

void Word::ToUpperCase( string & word )
{
    std::transform( word.begin(), word.end(), word.begin(), StringToUpperCase() );
}

EndNameSpace
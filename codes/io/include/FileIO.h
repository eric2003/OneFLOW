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
#include <fstream>
#include <sstream>
#include "Word.h"
#include "HXDefine.h"

BeginNameSpace( ONEFLOW )

const std::string whiteSpace = " ";

void SetDefaultLine( std::string * defaultLineIn );
std::string * GetDefaultLine();

void SetDefaultSeparatorOfWord( std::string * separatorOfWordIn );
std::string * GetDefaultSeparatorOfWord();
typedef std::streamsize StreamSize;

class CommentLine;
class FileIO
{
public:
    FileIO();
    ~FileIO();
protected:
    std::string * line, * separator;
    std::fstream file;
    int setfileFlag;
protected:
    std::string fileName;
    std::ios_base::openmode fileOpenMode;
    StreamSize filePosition;
    CommentLine * commentLine;
public:
    void OpenFile( const std::string & fileName, const std::ios_base::openmode & fileOpenMode );
    void OpenPrjFile( const std::string & fileName, const std::ios_base::openmode & fileOpenMode );
    void CloseFile();
    void MarkCurrentFilePosition();
    void MoveToPreviousFilePosition();
    std::string & GetCurrentLine() { return * this->line; };
public:
    void ResetCommentString( StringField &commentStringList );
    void SetDefaultSeparator( const std::string & separatorIn   ) { * this->separator = separatorIn; };

    std::string  * GetDefaultLine     () { return line; }
    std::fstream * GetDefaultFile     () { return &file; }
    std::string  * GetDefaultSeparator() { return separator; }
public:
    void SetLineContent( const std::string & lineContent ) { * this->line = lineContent; };
    void ShiftLineContent( int numberOfChars ) { * this->line = ( * this->line ).substr( numberOfChars ); }

    bool ReadNextMeaningfulLine();
    bool ReachTheEndOfFile();
public:
    void SkipLines( int numberOfLinesToSkip );
    bool ReadNextNonEmptyLine();
    bool NextWordIsEmpty();
    void DumpLineContentToScreen();
    std::string ReadNextWord();
    std::string ReadNextWord( const std::string & separator );
    std::string ReadNextTrueWord();
    std::string ReadNextWordToLowerCase();
    std::string ReadNextWordToLowerCase( const std::string & separator );
public:
    void SkipReadSymbol( const std::string & stringSymbol );
    void SkipReadWholeBlock();
public:
    template < typename T >
    friend inline FileIO & operator >> ( FileIO & textFileRead, T & value )
    {
        std::fstream & file = * textFileRead.GetDefaultFile();
        file >> value;
        return textFileRead;
    }

    template < typename T >
    T ReadNextDigit( std::ios_base & ( * f )( std::ios_base & ) = & std::dec )
    {
        std::string word = FileIO::ReadNextTrueWord();
        T value = StringToDigit< T >( word, f );
        return value;
    }

    template < typename T >
    T ReadNextDigit( int & num, std::ios_base & ( * f )( std::ios_base & ) = & std::dec )
    {
        std::string word = FileIO::ReadNextTrueWord();
        num = 1;
        bool flag = Word::FindString( word, "*" );
        if ( flag )
        {
            std::string word_num = Word::FindNextWord( word, "*" );
            num = StringToDigit< int >( word_num, f );
            word = Word::FindNextWord( word, "*" );
        }
        T value = StringToDigit< T >( word, f );
        return value;
    }

};

template < typename T >
inline T ReadNextDigit( std::ios_base & ( * f )( std::ios_base & ) = & std::dec )
{
    std::string * separatorOfWord = ONEFLOW::GetDefaultSeparatorOfWord();
    std::string * defaultLine     = ONEFLOW::GetDefaultLine();

    std::string word = Word::FindNextWord( * defaultLine, * separatorOfWord );
    T value = ONEFLOW::StringToDigit< T >( word, f );
    return value;
}

template < typename T >
inline T ReadNextDigit( std::string & source, const std::string & separatorOfWord, std::ios_base & ( * f )( std::ios_base & ) = & std::dec )
{
    std::string word = Word::FindNextWord( source, separatorOfWord );
    T value = ONEFLOW::StringToDigit< T >( word, f );
    return value;
}

bool IsEmpty( std::fstream & file );

EndNameSpace

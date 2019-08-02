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

#include "FileIO.h"
#include "CommentLine.h"
#include "Prj.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

std::string * defaultLine = 0;
void SetDefaultLine( std::string * defaultLineIn )
{
    defaultLine = defaultLineIn;
}

std::string * GetDefaultLine()
{
    return defaultLine;
}

std::string * separatorOfWord = 0;
void SetDefaultSeparatorOfWord( std::string * separatorOfWordIn )
{
    separatorOfWord = separatorOfWordIn;
}

std::string * GetDefaultSeparatorOfWord()
{
    return separatorOfWord;
}

FileIO::FileIO()
{
    line        = new std::string;
    separator   = new std::string;
    file        = new fstream;
    setfileFlag = 0;
    //\t is tab key
    string keyWordSeparator = " =\r\n\t#$,;\"";
    this->SetDefaultSeparator( keyWordSeparator );

    this->commentLine = new CommentLine();
    this->commentLine->AddString("#");
    this->commentLine->AddString("//");
}

FileIO::~FileIO()
{
    delete line;
    delete separator;
    if ( setfileFlag == 0 )
    {
        delete file;
    }
    delete this->commentLine;
}

void FileIO::ResetCommentString(StringField& commentStringList)
{
    this->commentLine->ResetCommentString(commentStringList);
}

void FileIO::SetDefaultFile ( std::fstream * defaultFileIn )
{
    if ( setfileFlag == 0 )
    {
        delete file;
    }
    this->file  = defaultFileIn;
    setfileFlag = 1;
}

void FileIO::OpenFile( const string & fileName, const ios_base::openmode & fileOpenMode )
{
    this->fileName     = fileName;
    this->fileOpenMode = fileOpenMode;
    ONEFLOW::OpenFile( * file, fileName, fileOpenMode );
}

void FileIO::OpenPrjFile( const string & fileName, const ios_base::openmode & fileOpenMode )
{
    this->fileName     = fileName;
    this->fileOpenMode = fileOpenMode;
    ONEFLOW::OpenPrjFile( * file, fileName, fileOpenMode );
}

void FileIO::CloseFile()
{
    ONEFLOW::CloseFile( * file );
}

void FileIO::MarkCurrentFilePosition()
{
    filePosition = file->tellp();
}

void FileIO::MoveToPreviousFilePosition()
{
    file->seekp( filePosition );
}

bool FileIO::ReadNextMeaningfulLine()
{
    while ( ! this->ReachTheEndOfFile() )
    {
        ONEFLOW::ReadNextLine( * file, * line );

        if ( ONEFLOW::IsEmptyLine  ( * line ) ||
             ONEFLOW::IsCommentLine( * line, this->commentLine->commentdata ) )
        {
             continue;
        }
        return true;
    }
    return false;
}

bool FileIO::ReachTheEndOfFile()
{
    if ( ( * file ).eof() )
    {
        return true;
    }
    return false;
}

void FileIO::SkipLines( int numberOfLinesToSkip )
{
    ONEFLOW::SkipLines( * file, numberOfLinesToSkip );
}

bool FileIO::ReadNextNonEmptyLine()
{
    return ONEFLOW::ReadNextNonEmptyLine( * this->file, * this->line );
}

void FileIO::DumpLineContentToScreen()
{
    cout << * line << endl;
}

void FileIO::SkipReadSymbol( const string & stringSymbol )
{
    while ( ! this->ReachTheEndOfFile() )
    {
        bool resultFlag = this->ReadNextMeaningfulLine();
        if ( ! resultFlag ) break;

        string word = this->ReadNextWord();

        if ( word == stringSymbol )
        {
            return;
        }
    }
}

void FileIO::SkipReadWholeBlock()
{
    int countOfLeftBrackets  = 0;
    int countOfRightBrackets = 0;

    while ( ! this->ReachTheEndOfFile() )
    {
        bool resultFlag = this->ReadNextMeaningfulLine();
        if ( ! resultFlag ) break;

        string word = this->ReadNextWord();

        if ( word == "{" )
        {
            ++ countOfLeftBrackets;
        }
        else if ( word == "}" )
        {
            ++ countOfRightBrackets;
        }

        if ( countOfLeftBrackets == countOfRightBrackets )
        {
            return;
        }
    }
}

bool FileIO::NextWordIsEmpty()
{
    std::string lineLeft = * this->line;
    string word = ONEFLOW::FindNextWord( lineLeft, * this->separator );
    return word == "";
}

std::string FileIO::ReadNextTrueWord()
{
    string word = ONEFLOW::FindNextWord( * this->line, * this->separator );

    if ( word == "" )
    {
        ReadNextNonEmptyLine();
        word = ONEFLOW::FindNextWord( * this->line, * this->separator );
    }

    return word;
}

std::string FileIO::ReadNextWord()
{
    string word = ONEFLOW::FindNextWord( * this->line, * this->separator );

    return word;
}

std::string FileIO::ReadNextWord( const std::string & separator )
{
    string word = ONEFLOW::FindNextWord( * this->line, separator );
    return word;
}

string FileIO::ReadNextWordToLowerCase()
{
    string word = ONEFLOW::FindNextWord( * this->line, * this->separator );
    ONEFLOW::ToLowerCase( word );
    return word;
}

string FileIO::ReadNextWordToLowerCase( const std::string & separator )
{
    string word = ONEFLOW::FindNextWord( * this->line, separator );
    ONEFLOW::ToLowerCase( word );
    return word;
}

bool IsEmpty( fstream & file )
{
    file.seekp( 0, ios::end );
    streamoff i = file.tellp();
    //cout << "ONEFLOW::IsEmpty( fstream & file ) = " << i << endl;
    if ( i ) return false;
    return true;
}

bool IsEmptyLine( const string & line )
{
    if ( line == "" )
    {
        return true;
    }
    else
    {
        const string notReadableSeparator = " \r\n\t";
        string word;
        FindNextWord( line, word, notReadableSeparator );
        return word == "";
    }
}

bool IsCommentLine( const string & line )
{
    const string notReadableSeparator = " \r\n\t";
    string word;
    FindNextWord( line, word, notReadableSeparator );
    return ( word.substr( 0, 1 ) == "#" ||
        word.substr( 0, 2 ) == "//" );
}

bool IsCommentLine(const string& line, StringField &comlist )
{
    const string notReadableSeparator = " \r\n\t";
    string word;
    FindNextWord(line, word, notReadableSeparator);
    for (int i = 0; i < comlist.size(); ++ i)
    {
        string & t = comlist[ i ];
        int n = t.length();
        if ( word.substr(0, n) == t ) return true;
    }
    return false;
}

EndNameSpace
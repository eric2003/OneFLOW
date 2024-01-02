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

#include "FileIO.h"

#include "Word.h"
#include "CommentLine.h"
#include "Prj.h"
#include <iostream>

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
    line        = new std::string();
    separator   = new std::string();
    //file        = new std::fstream();
    setfileFlag = 0;
    //\t is tab key
    std::string keyWordSeparator = " =\r\n\t#$,;\"";
    this->SetDefaultSeparator( keyWordSeparator );

    this->commentLine = new CommentLine();
    this->commentLine->AddString("#");
    this->commentLine->AddString("//");
}

FileIO::~FileIO()
{
    delete line;
    delete separator;
    delete this->commentLine;
}

void FileIO::ResetCommentString(StringField& commentStringList)
{
    this->commentLine->ResetCommentString(commentStringList);
}

void FileIO::OpenFile( const std::string & fileName, const std::ios_base::openmode & fileOpenMode )
{
    this->fileName     = fileName;
    this->fileOpenMode = fileOpenMode;
    Prj::OpenFile( this->file, fileName, fileOpenMode );
}

void FileIO::OpenPrjFile( const std::string & fileName, const std::ios_base::openmode & fileOpenMode )
{
    this->fileName     = fileName;
    this->fileOpenMode = fileOpenMode;
    Prj::OpenPrjFile( this->file, fileName, fileOpenMode );
}

void FileIO::CloseFile()
{
    Prj::CloseFile( this->file );
}

void FileIO::MarkCurrentFilePosition()
{
    filePosition = this->file.tellp();
}

void FileIO::MoveToPreviousFilePosition()
{
    this->file.seekp( filePosition );
}

bool FileIO::ReadNextMeaningfulLine()
{
    while ( ! this->ReachTheEndOfFile() )
    {
         Word::ReadNextLine( this->file, * line );

        if ( Word::IsEmptyLine  ( * line ) ||
             Word::IsCommentLine( * line, this->commentLine->commentdata ) )
        {
             continue;
        }
        return true;
    }
    return false;
}

bool FileIO::ReachTheEndOfFile()
{
    if ( this->file.eof() )
    {
        return true;
    }
    return false;
}

void FileIO::SkipLines( int numberOfLinesToSkip )
{
     Word::SkipLines( this->file, numberOfLinesToSkip );
}

bool FileIO::ReadNextNonEmptyLine()
{
    return Word::ReadNextNonEmptyLine( this->file, * this->line );
}

void FileIO::DumpLineContentToScreen()
{
    std::cout << * line << std::endl;
}

void FileIO::SkipReadSymbol( const std::string & stringSymbol )
{
    while ( ! this->ReachTheEndOfFile() )
    {
        bool resultFlag = this->ReadNextMeaningfulLine();
        if ( ! resultFlag ) break;

        std::string word = this->ReadNextWord();

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

        std::string word = this->ReadNextWord();

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
    std::string word = Word::FindNextWord( lineLeft, * this->separator );
    return word == "";
}

std::string FileIO::ReadNextTrueWord()
{
    std::string word = Word::FindNextWord( * this->line, * this->separator );

    if ( word == "" )
    {
        ReadNextNonEmptyLine();
        word = Word::FindNextWord( * this->line, * this->separator );
    }

    return word;
}

std::string FileIO::ReadNextWord()
{
    std::string word = Word::FindNextWord( * this->line, * this->separator );

    return word;
}

std::string FileIO::ReadNextWord( const std::string & separator )
{
    std::string word = Word::FindNextWord( * this->line, separator );
    return word;
}

std::string FileIO::ReadNextWordToLowerCase()
{
    std::string word = Word::FindNextWord( * this->line, * this->separator );
    Word::ToLowerCase( word );
    return word;
}

std::string FileIO::ReadNextWordToLowerCase( const std::string & separator )
{
    std::string word = Word::FindNextWord( * this->line, separator );
    Word::ToLowerCase( word );
    return word;
}

bool IsEmpty( std::fstream & file )
{
    file.seekp( 0, std::ios::end );
    std::streamoff i = file.tellp();
    //std::cout << "ONEFLOW::IsEmpty( std::fstream & file ) = " << i << std::endl;
    if ( i ) return false;
    return true;
}

EndNameSpace

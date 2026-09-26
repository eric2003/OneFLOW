/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2026 He Xin and the OneFLOW contributors.
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

#include "TextFileParser.h"

#include "Word.h"
#include "Prj.h"
#include <iostream>

BeginNameSpace( ONEFLOW )

TextFileParser::TextFileParser()
{
    // \t is a literal tab character.
    separator = " =\r\n\t#$,;\"";

    commentLine.AddString( "#" );
    commentLine.AddString( "//" );
}

void TextFileParser::ResetCommentString( StringField & commentStringList )
{
    commentLine.ResetCommentString( commentStringList );
}

void TextFileParser::OpenFile( const std::string & fileName, const std::ios_base::openmode & fileOpenMode )
{
    Prj::OpenFile( file, fileName, fileOpenMode );
}

void TextFileParser::OpenPrjFile( const std::string & fileName, const std::ios_base::openmode & fileOpenMode )
{
    Prj::OpenPrjFile( file, fileName, fileOpenMode );
}

void TextFileParser::CloseFile()
{
    Prj::CloseFile( file );
}

void TextFileParser::MarkCurrentFilePosition()
{
    filePosition = file.tellp();
}

void TextFileParser::MoveToPreviousFilePosition()
{
    file.seekp( filePosition );
}

bool TextFileParser::ReadNextMeaningfulLine()
{
    while ( ! ReachTheEndOfFile() )
    {
        Word::ReadNextLine( file, line );

        if ( Word::IsEmptyLine( line ) ||
            Word::IsCommentLine( line, commentLine.commentdata ) )
        {
            continue;
        }
        return true;
    }
    return false;
}

bool TextFileParser::ReachTheEndOfFile()
{
    return file.eof();
}

void TextFileParser::SkipLines( int numberOfLinesToSkip )
{
    Word::SkipLines( file, numberOfLinesToSkip );
}

bool TextFileParser::ReadNextNonEmptyLine()
{
    return Word::ReadNextNonEmptyLine( file, line );
}

void TextFileParser::DumpLineContentToScreen()
{
    std::cout << line << std::endl;
}

void TextFileParser::SkipReadSymbol( const std::string & stringSymbol )
{
    while ( ! ReachTheEndOfFile() )
    {
        if ( ! ReadNextMeaningfulLine() ) break;

        if ( ReadNextWord() == stringSymbol )
        {
            return;
        }
    }
}

void TextFileParser::SkipReadWholeBlock()
{
    int countOfLeftBrackets  = 0;
    int countOfRightBrackets = 0;

    while ( ! ReachTheEndOfFile() )
    {
        if ( ! ReadNextMeaningfulLine() ) break;

        std::string word = ReadNextWord();

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

bool TextFileParser::NextWordIsEmpty()
{
    std::string lineLeft = line;
    return Word::FindNextWord( lineLeft, separator ).empty();
}

std::string TextFileParser::ReadNextTrueWord()
{
    std::string word = Word::FindNextWord( line, separator );

    if ( word.empty() )
    {
        ReadNextNonEmptyLine();
        word = Word::FindNextWord( line, separator );
    }

    return word;
}

std::string TextFileParser::ReadNextWord()
{
    return Word::FindNextWord( line, separator );
}

std::string TextFileParser::ReadNextWord( const std::string & separatorIn )
{
    return Word::FindNextWord( line, separatorIn );
}

std::string TextFileParser::ReadNextWordToLowerCase()
{
    std::string word = Word::FindNextWord( line, separator );
    Word::ToLowerCase( word );
    return word;
}

std::string TextFileParser::ReadNextWordToLowerCase( const std::string & separatorIn )
{
    std::string word = Word::FindNextWord( line, separatorIn );
    Word::ToLowerCase( word );
    return word;
}

bool IsEmpty( std::fstream & file )
{
    file.seekp( 0, std::ios::end );
    return file.tellp() == 0;
}

EndNameSpace


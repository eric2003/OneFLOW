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
#pragma once

#include <fstream>
#include <ios>
#include <sstream>
#include <string>

#include "Word.h"
#include "HXDefine.h"
#include "CommentLine.h"

BeginNameSpace( ONEFLOW )

// --- Dead code removed in this refactor (verified against all ~810 files) ---
// * Free functions SetDefaultLine/GetDefaultLine/SetDefaultSeparatorOfWord/
//   GetDefaultSeparatorOfWord, and the zero-argument free ReadNextDigit<T>()
//   templates that read from them: these mutated file-scope global pointers
//   and had ZERO call sites anywhere outside their own declaration/definition.
// * Member accessors GetDefaultLine()/GetDefaultFile()/GetDefaultSeparator():
//   zero call sites (GetDefaultFile() was only used by the dead operator>>
//   below).
// * SetLineContent(), ShiftLineContent(): declared, never called.
// * The friend template operator>>: declared, never called.
// * The `whiteSpace` constant: declared, never referenced.
// * The `setfileFlag` data member: set to 0 in the constructor, never read.
// * The `fileName` / `fileOpenMode` data members: assigned in OpenFile()/
//   OpenPrjFile(), never read back by anything (Prj::OpenFile/OpenPrjFile
//   receive the filename/mode as ordinary parameters, not through `this`).
// Removing these deletes two raw global pointers and ~40 lines of code that
// no longer had any callers, without changing behavior for any real caller.

class TextFileParser
{
public:
    TextFileParser();
    ~TextFileParser() = default;

    // std::fstream is not copyable, so this type was already implicitly
    // non-copyable; that is now explicit. Move is intentionally not enabled:
    // nothing in the codebase moves a TextFileParser, it is always a local,
    // stack-lived object opened/read/closed in one place.
    TextFileParser( const TextFileParser & ) = delete;
    TextFileParser & operator=( const TextFileParser & ) = delete;

public:
    void OpenFile( const std::string & fileName, const std::ios_base::openmode & fileOpenMode );
    void OpenPrjFile( const std::string & fileName, const std::ios_base::openmode & fileOpenMode );
    void CloseFile();

    void MarkCurrentFilePosition();
    void MoveToPreviousFilePosition();

    [[nodiscard]] std::string & GetCurrentLine() { return line; }

public:
    void ResetCommentString( StringField & commentStringList );
    void SetDefaultSeparator( const std::string & separatorIn ) { separator = separatorIn; }

    [[nodiscard]] bool ReadNextMeaningfulLine();
    [[nodiscard]] bool ReachTheEndOfFile();

public:
    void SkipLines( int numberOfLinesToSkip );
    bool ReadNextNonEmptyLine();
    [[nodiscard]] bool NextWordIsEmpty();
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
    T ReadNextDigit( std::ios_base & ( * f )( std::ios_base & ) = & std::dec )
    {
        std::string word = ReadNextTrueWord();
        return StringToDigit< T >( word, f );
    }

    // Handles the "N*value" repeat-count shorthand used by several
    // system/*.txt config files (e.g. "5*0.0" == the value 0.0, five times).
    // `num` receives the repeat count (1 if the shorthand wasn't used).
    template < typename T >
    T ReadNextDigit( int & num, std::ios_base & ( * f )( std::ios_base & ) = & std::dec )
    {
        std::string word = ReadNextTrueWord();
        num = 1;

        if ( Word::FindString( word, "*" ) )
        {
            std::string repeatCountWord = Word::FindNextWord( word, "*" );
            num  = StringToDigit< int >( repeatCountWord, f );
            word = Word::FindNextWord( word, "*" );
        }

        return StringToDigit< T >( word, f );
    }

private:
    std::string line;
    std::string separator;
    std::fstream file;
    CommentLine commentLine;

    // Only meaningful between a MarkCurrentFilePosition() /
    // MoveToPreviousFilePosition() pair.
    std::streamsize filePosition = 0;
};

[[nodiscard]] bool IsEmpty( std::fstream & file );

EndNameSpace


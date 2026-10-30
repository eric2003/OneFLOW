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

#include "ParaFile.h"
#include "TextFileParser.h"
#include "Parallel.h"
#include "Prj.h"
#include <ios>
#include <string>
#include <vector>

BeginNameSpace( ONEFLOW )

void ReadControlInfo()
{
    if ( Parallel::IsServer() )
    {
        ONEFLOW::ReadPrjScript();
    }

    Parallel::TestSayHelloFromEveryProcess();
    ONEFLOW::BroadcastControlParameterToAllProcessors();
    ONEFLOW::DumpDataBase();
}

void ReadControlInfo( const std::string & caseDir )
{
    if ( Parallel::IsServer() )
    {
        ONEFLOW::ReadPrjScript( caseDir );
    }

    Parallel::TestSayHelloFromEveryProcess();
    ONEFLOW::BroadcastControlParameterToAllProcessors();
    ONEFLOW::DumpDataBase( caseDir );
}

void ReadPrjScript()
{
    std::vector< std::string > scriptFileNameList;
    ONEFLOW::ReadScriptFileNameList( scriptFileNameList );
    ONEFLOW::ReadMultiScriptFiles( scriptFileNameList );
}

void ReadPrjScript( const std::string & caseDir )
{
    std::vector< std::string > scriptFileNameList;
    ONEFLOW::ReadScriptFileNameList( caseDir, scriptFileNameList );
    ONEFLOW::ReadMultiScriptFiles( scriptFileNameList );
}

void ReadScriptFileNameList( std::vector< std::string > & scriptFileNameList )
{
    TextFileParser textFileParser;

    textFileParser.OpenFile(
        Prj::GetPrjFileName( "script/control.txt" ),
        std::ios_base::in );

    // Tab is a separator.
    std::string keyWordSeparator = " ()\r\n\t#$,;\"";
    textFileParser.SetDefaultSeparator( keyWordSeparator );

    while ( ! textFileParser.ReachTheEndOfFile() )
    {
        bool flag = textFileParser.ReadNextNonEmptyLine();
        if ( ! flag ) break;

        std::string scriptFileName =
            textFileParser.ReadNextWord();

        std::string fullScriptFileName =
            Prj::GetPrjFileName(
                "script/" + scriptFileName );

        scriptFileNameList.push_back( fullScriptFileName );
    }

    textFileParser.CloseFile();
}

void ReadScriptFileNameList(
    const std::string & caseDir,
    std::vector< std::string > & scriptFileNameList )
{
    TextFileParser textFileParser;
    textFileParser.OpenCaseFile(
        caseDir,
        "script/control.txt",
        std::ios_base::in );

    // Tab is a separator.
    std::string keyWordSeparator = " ()\r\n\t#$,;\"";
    textFileParser.SetDefaultSeparator( keyWordSeparator );

    while ( ! textFileParser.ReachTheEndOfFile() )
    {
        bool flag = textFileParser.ReadNextNonEmptyLine();
        if ( ! flag ) break;

        std::string scriptFileName =
            textFileParser.ReadNextWord();

        std::string fullScriptFileName =
            Prj::GetCaseFileName(
                caseDir,
                "script/" + scriptFileName );

        scriptFileNameList.push_back( fullScriptFileName );
    }

    textFileParser.CloseFile();
}

void ReadMultiScriptFiles( std::vector< std::string > & scriptFileNameList )
{
    for ( const std::string & scriptFileName : scriptFileNameList )
    {
        ONEFLOW::ReadOneFLOWScriptFile( scriptFileName );
    }
}

EndNameSpace

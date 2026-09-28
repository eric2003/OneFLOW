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

#include "SolverNameList.h"
#include "SolverNamePolicy.h"
#include "TextFileParser.h"  // only if needed in header; prefer .cpp only
#include "GridState.h"

BeginNameSpace( ONEFLOW )


StringField SolverNameClass::unsSolverNameList;
StringField SolverNameClass::strSolverNameList;
bool SolverNameClass::flag = false;

SolverNameClass::SolverNameClass()
{
    ;
}

SolverNameClass::~SolverNameClass()
{
    ;
}

void SolverNameClass::Init()
{
    if ( flag ) return;
    flag = true;
    SolverNameClass::ReadSolverNames();
}

void SolverNameClass::ReadSolverNames()
{
    StringField solverNameList;
    SolverNameClass::ReadSolverNames( solverNameList );

    FillExpandedSolverNames(
        solverNameList,
        SolverNameClass::unsSolverNameList,
        SolverNameClass::strSolverNameList );
}

void SolverNameClass::ReadSolverNames( StringField & solverNameList )
{
    TextFileParser textFileParser;

    textFileParser.OpenPrjFile( "script/solver.txt", std::ios_base::in );

    // \t is the tab key
    const std::string keyWordSeparator = " ()\r\n\t#$,;\"";
    textFileParser.SetDefaultSeparator( keyWordSeparator );

    // Same pattern as MessageMapImp::ReadFile:
    // skip blank/comment lines; no spurious empty token at EOF.
    while ( textFileParser.ReadNextMeaningfulLine() )
    {
        std::string solverName = textFileParser.ReadNextWord();
        if ( solverName.empty() )
        {
            continue;
        }
        solverNameList.push_back( solverName );
    }

    textFileParser.CloseFile();
}

StringField & SolverNameClass::GetSolverNames( int gridType )
{
    SolverNameClass::Init();

    if ( gridType == ONEFLOW::UMESH )
    {
        return SolverNameClass::unsSolverNameList;
    }
    else
    {
        return SolverNameClass::strSolverNameList;
    }
}

void SolverNameClass::LoadFromBaseNames( const StringField & baseNames )
{
    SolverNameClass::unsSolverNameList.clear();
    SolverNameClass::strSolverNameList.clear();

    FillExpandedSolverNames(
        baseNames,
        SolverNameClass::unsSolverNameList,
        SolverNameClass::strSolverNameList );

    flag = true; // skip file path on later Init()
}

void SolverNameClass::Reset()
{
    unsSolverNameList.clear();
    strSolverNameList.clear();
    flag = false;
}

EndNameSpace
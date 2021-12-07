/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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
#include "FileIO.h"
#include "DataBase.h"
#include "DataBook.h"
#include <vector>
#include <string>
using namespace std;

BeginNameSpace( ONEFLOW )

bool IsArrayParameter( const std::string & lineOfName );
void ReadOneFLOWScriptFile( FileIO & fileIO );
void ReadOneFLOWScriptFile( const std::string & fileName );
string GetJsonFileName( const std::string & fileName );
void GetParaInfo( FileIO & fileIO, std::string & varName, std::vector< std::string > & varArray );
void GetParaInfoArray( FileIO & fileIO, std::string & varName, std::vector< std::string > & varArray );
void GetParaInfoScalar( FileIO & fileIO, std::string & varName, std::vector< std::string > & varArray );

void AnalysisArrayParameter( FileIO & fileIO, int keyWordIndex );
int AnalysisScalarParameter( FileIO & fileIO, int keyWordIndex );
int GetParameterArraySize( const std::string & word );

void ReadControlInfo();
void ReadPrjScript();
void ReadScriptFileNameList( std::vector< std::string > & scriptFileNameList );
void ReadMultiScriptFiles( std::vector< std::string > & scriptFileNameList );
void BroadcastControlParameterToAllProcessors();
void DumpDataBase();

void CompressData( DataBase * dataBase, DataBook *& dataBook );
void DecompressData( DataBase * dataBase, DataBook * dataBook );

void CompressData( DataBook *& dataBook );
void DecompressData( DataBook * dataBook );

EndNameSpace

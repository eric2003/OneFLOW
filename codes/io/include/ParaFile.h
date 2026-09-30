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
#include <vector>
#include <string>

namespace ONEFLOW {

class DataBase;
class DataBook;
class TextFileParser;

bool IsArrayParameter( const std::string & lineOfName );
void ReadOneFLOWScriptFile( TextFileParser & textFileParser );
void ReadOneFLOWScriptFile( const std::string & fileName );
std::string GetJsonFileName( const std::string & fileName );
void GetParaInfo( TextFileParser & textFileParser, std::string & varName, std::vector< std::string > & varArray );
void GetParaInfoArray( TextFileParser & textFileParser, std::string & varName, std::vector< std::string > & varArray );
void GetParaInfoScalar( TextFileParser & textFileParser, std::string & varName, std::vector< std::string > & varArray );

void AnalysisArrayParameter( TextFileParser & textFileParser, int keyWordIndex );
int AnalysisScalarParameter( TextFileParser & textFileParser, int keyWordIndex );
int GetParameterArraySize( const std::string & word );

void ReadControlInfo();
void ReadControlInfo( const std::string & caseDir );
void ReadPrjScript();
void ReadPrjScript( const std::string & caseDir );
void ReadScriptFileNameList( std::vector< std::string > & scriptFileNameList );
void ReadScriptFileNameList( const std::string & caseDir, std::vector< std::string > & scriptFileNameList );
void ReadMultiScriptFiles( std::vector< std::string > & scriptFileNameList );
void BroadcastControlParameterToAllProcessors();
void DumpDataBase();
void DumpDataBase( const std::string & caseDir );

void CompressData( DataBase * dataBase, DataBook *& dataBook );
void DecompressData( DataBase * dataBase, DataBook * dataBook );

void CompressData( DataBook *& dataBook );
void DecompressData( DataBook * dataBook );

} // namespace ONEFLOW

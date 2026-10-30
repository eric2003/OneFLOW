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

/**
 * @brief Loads project control scripts and completes process-level setup.
 *
 * ReadControlInfo owns the parallel greeting, control-parameter broadcast,
 * and database diagnostic output in addition to server-side script loading.
 * Parallel state must be initialized before this function is called.
 */
void ReadControlInfo();
void ReadControlInfo( const std::string & caseDir );

/**
 * @brief Discovers and commits project control scripts in manifest order.
 *
 * These functions perform script ingestion only. They do not initialize
 * parallel state, broadcast parameters, or write the runtime database dump.
 */
void ReadPrjScript();
void ReadPrjScript( const std::string & caseDir );

/**
 * @brief Appends script paths from the ordered control manifest.
 *
 * The case overload resolves entries under the supplied case directory.
 * Existing vector entries are preserved; manifest order is preserved.
 */
void ReadScriptFileNameList( std::vector< std::string > & scriptFileNameList );
void ReadScriptFileNameList( const std::string & caseDir, std::vector< std::string > & scriptFileNameList );
/**
 * @brief Reads script files in manifest order and commits each document.
 *
 * Each file is parsed and committed before the next file is read. The
 * operation is ordered but not transactional across the complete file list.
 */
void ReadMultiScriptFiles( const std::vector< std::string > & scriptFileNameList );
void BroadcastControlParameterToAllProcessors();
void DumpDataBase();
void DumpDataBase( const std::string & caseDir );

void CompressData( DataBase * dataBase, DataBook * dataBook );
void DecompressData( DataBase * dataBase, DataBook * dataBook );

void CompressData( DataBook * dataBook );
void DecompressData( DataBook * dataBook );

} // namespace ONEFLOW

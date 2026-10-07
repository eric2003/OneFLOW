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
#include <memory>
#include "TextFileParser.h"
#include "DataBase.h"
#include "DataBook.h"
#include "ConfigLoader.h"
#include "ConfigDatabaseAdapter.h"
#include "LegacyParameterSyntax.h"
#include "DataBase.h"
#include "Parallel.h"
#include "LogFile.h"
#include "OStream.h"
#include "Fatal.h"
#include "Prj.h"
#include "FileUtils.h"
#include "PIO.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>


BeginNameSpace( ONEFLOW )

bool IsArrayParameter( const std::string & lineOfName )
{
    return ONEFLOW::IsLegacyArrayParameter( lineOfName );
}

int GetParameterArraySize( const std::string & word )
{
    if ( Word::IsDigit( word ) )
    {
        return StringToDigit< int >( word );
    }
    return GetDataValue< int >( word );
}

void ReadOneFLOWScriptFile( TextFileParser & textFileParser )
{
    ConfigLoader loader( GetParameterArraySize );
    loader.ParseFromParser( textFileParser );
    ConfigDatabaseAdapter::Commit( loader.Document() );
}

void ReadOneFLOWScriptFile( const std::string & fileName )
{
    ConfigLoader loader( GetParameterArraySize );
    loader.ParseFile( fileName );
    ConfigDatabaseAdapter::Commit( loader.Document() );
}

void AnalysisArrayParameter( TextFileParser & textFileParser, int keyWordIndex )
{
    std::string errorMessage = "error in parameter file";
    std::string commSeparator = "=\r\n\t#$,;\"";

    std::string ayrrayInfo = textFileParser.ReadNextWord( commSeparator );

    //Array pattern
    std::string arraySeparator = " =\r\n\t#$,;\"[]";
    std::string arrayName, arraySizeName;

    arrayName = Word::FindNextWord( ayrrayInfo, arraySeparator );
    arraySizeName = Word::FindNextWord( ayrrayInfo, arraySeparator );

    int arraySize = ONEFLOW::GetParameterArraySize( arraySizeName );

    std::vector<std::string> valueContainer( static_cast<std::size_t>( arraySize ) );

    for ( int i = 0; i < arraySize; ++ i )
    {
        valueContainer[ i ] = textFileParser.ReadNextWord( arraySeparator );
        //It shows that these contents can't be written within 1 lines
        if ( valueContainer[ i ] == "" )
        {
            textFileParser.ReadNextNonEmptyLine();
            valueContainer[ i ] = textFileParser.ReadNextWord( arraySeparator );
            if ( valueContainer[ i ] == "" )
            {
                Fatal( errorMessage );
            }
        }
    }
    ONEFLOW::ProcessData( arrayName, valueContainer.data(), keyWordIndex, arraySize );

}

int AnalysisScalarParameter( TextFileParser & textFileParser, int keyWordIndex )
{
    std::string errorMessage = "error in parameter file";
    std::string separator = " =\r\n\t#$,;\"";  //\t is tab key

    std::string name = textFileParser.ReadNextWord( separator );

    int arraySize = 1;
    std::vector<std::string> value( static_cast<std::size_t>( arraySize ) );

    value[ 0 ] = textFileParser.ReadNextWord( separator );

    ONEFLOW::ProcessData( name, value.data(), keyWordIndex, arraySize );

    return arraySize;
}

std::string GetJsonFileName( const std::string & fileName )
{
    std::string mainName, extensionName;
    ONEFLOW::GetFileNameExtension( fileName, mainName, extensionName, "." );
    std::string newExtensionName = "json";

    OStream &logger = OStream::Instance();
    logger.ClearAll();
    logger << mainName << "." << newExtensionName;

    std::string newFileName = logger.str();
    return newFileName;
}

void GetParaInfo( TextFileParser & textFileParser, std::string & varName, std::vector< std::string > & varArray )
{
    if ( ONEFLOW::IsArrayParameter( textFileParser.GetCurrentLine() ) )
    {
        ONEFLOW::GetParaInfoArray( textFileParser, varName, varArray );
    }
    else
    {
        ONEFLOW::GetParaInfoScalar( textFileParser, varName, varArray );
    }
}

void GetParaInfoScalar( TextFileParser & textFileParser, std::string & varName, std::vector< std::string > & varArray )
{
    std::string errorMessage = "error in parameter file";
    std::string separator = " =\r\n\t#$,;\"";  //\t is tab key

    varName = textFileParser.ReadNextWord( separator );

    int arraySize = 1;
    varArray.resize( arraySize );

    varArray[ 0 ] = textFileParser.ReadNextWord( separator );
}

void GetParaInfoArray( TextFileParser & textFileParser, std::string & varName, std::vector< std::string > & varArray )
{
    std::string errorMessage = "error in parameter file";
    std::string commSeparator = "=\r\n\t#$,;\"";

    std::string ayrrayInfo = textFileParser.ReadNextWord( commSeparator );

    //Array pattern
    std::string arraySeparator = " =\r\n\t#$,;\"[]";
    std::string arrayName, arraySizeName;

    arrayName = Word::FindNextWord( ayrrayInfo, arraySeparator );
    arraySizeName = Word::FindNextWord( ayrrayInfo, arraySeparator );

    int arraySize = ONEFLOW::GetParameterArraySize( arraySizeName );

    varArray.resize( arraySize );

    for ( int i = 0; i < arraySize; ++ i )
    {
        varArray[ i ] = textFileParser.ReadNextWord( arraySeparator );
        //It shows that these contents can't be written within 1 lines
        if ( varArray[ i ] == "" )
        {
            textFileParser.ReadNextNonEmptyLine();
            varArray[ i ] = textFileParser.ReadNextWord( arraySeparator );
            if ( varArray[ i ] == "" )
            {
                Fatal( errorMessage );
            }
        }
    }
}

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

void DumpDataBase()
{
    DataBase & dataBase = ONEFLOW::RequireGlobalDataBase();
    std::fstream file;
    std::string fileName = "/log/database.log";
    PIO::OpenPrjFile( file, fileName, std::ios_base::out );
    dataBase.RequireDataPara().DumpData( file );
    PIO::CloseFile( file );
}

void DumpDataBase( const std::string & caseDir )
{
    DataBase * dataBase = ONEFLOW::GetGlobalDataBase();
    std::fstream file;
    Prj::OpenCaseFile( file, caseDir, "log/database.log", std::ios_base::out );
    dataBase->GetDataPara()->DumpData( file );
    PIO::CloseFile( file );
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

    textFileParser.OpenPrjFile(
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
    int numberOfParameterFiles = scriptFileNameList.size();

    for ( int iFile = 0; iFile < numberOfParameterFiles; ++ iFile )
    {
        std::string & scriptFileName = scriptFileNameList[ iFile ];

        ONEFLOW::ReadOneFLOWScriptFile( scriptFileName );
    }
}

void BroadcastControlParameterToAllProcessors()
{
    ONEFLOW::logFile << "Broadcast Control Parameter To All Processors\n";

    ONEFLOW::HXBcast( ONEFLOW::CompressData, ONEFLOW::DecompressData, Parallel::GetServerid() );
}

void CompressData( DataBook * dataBook )
{
    DataBase & globalDataBase = ONEFLOW::RequireGlobalDataBase();

    ONEFLOW::CompressData( &globalDataBase, dataBook );
}

void DecompressData( DataBook * dataBook )
{
    DataBase & globalDataBase = ONEFLOW::RequireGlobalDataBase();
    ONEFLOW::DecompressData( &globalDataBase, dataBook );
}

void CompressData( DataBase * dataBase, DataBook * dataBook )
{
    // Use the new type alias
    const DataPara::DataMap & dataMap = dataBase->GetDataPara()->GetDataMap();

    int ndata = static_cast<int>( dataMap.size() );
    ONEFLOW::HXWrite( dataBook, ndata );

    // Range-based for is cleaner with unordered_map
    for ( const auto & pair : dataMap )
    {
        const DataEntry * dataEntry = pair.second.get();  // pair.first is the key (name), pair.second owns DataEntry
        ONEFLOW::HXWriteDataEntry( dataBook, dataEntry );
    }
}

void DecompressData( DataBase * dataBase, DataBook * dataBook )
{
    // No longer need to touch the internal map directly for reading
    dataBook->MoveToBegin();

    int ndata = 0;
    ONEFLOW::HXRead( dataBook, ndata );

    for ( int i = 0; i < ndata; ++ i )
    {
        auto dataEntry = ONEFLOW::HXReadDataEntry( dataBook );
        dataBase->GetDataPara()->SetDataEntry( std::move( dataEntry ) );
    }
}

EndNameSpace

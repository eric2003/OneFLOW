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
#include "DataBase.h"
#include "Parallel.h"
#include "LogFile.h"
#include "OStream.h"
#include "Stop.h"
#include "Prj.h"
#include "FileUtils.h"
#include "PIO.h"
#include "json/json.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>


BeginNameSpace( ONEFLOW )

bool IsArrayParameter( const std::string & lineOfName )
{
    const std::string::size_type npos = - 1;

    if ( lineOfName.find_first_of( "[" ) == npos )
    {
        return false;
    }

    if ( lineOfName.find_first_of( "]" ) == npos )
    {
        return false;
    }

    return true;
}

void ReadOneFLOWScriptFile( TextFileParser & textFileParser )
{
    //string name, word;

    //\t is the tab key
    std::string keyWordSeparator = " =\r\n\t#$,;\"";

    textFileParser.SetDefaultSeparator( keyWordSeparator );

    DataBaseType::Init();

    while ( ! textFileParser.ReachTheEndOfFile() )
    {
        bool resultFlag = textFileParser.ReadNextMeaningfulLine();
        if ( ! resultFlag ) break;

        std::string keyWord = textFileParser.ReadNextWord();

        if ( keyWord == "" ) continue;

        //int keyWordIndex = keyWordMap[ keyWord ];
        int keyWordIndex = DataBaseType::GetIndex( keyWord );

        if ( ONEFLOW::IsArrayParameter( textFileParser.GetCurrentLine() ) )
        {
            ONEFLOW::AnalysisArrayParameter( textFileParser, keyWordIndex );
        }
        else
        {
            ONEFLOW::AnalysisScalarParameter( textFileParser, keyWordIndex );
        }
    }
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

    std::string * valueContainer = new std::string[ arraySize ];

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
                Stop( errorMessage );
            }
        }
    }
    ONEFLOW::ProcessData( arrayName, valueContainer, keyWordIndex, arraySize );

    delete[] valueContainer;
}

int AnalysisScalarParameter( TextFileParser & textFileParser, int keyWordIndex )
{
    std::string errorMessage = "error in parameter file";
    std::string separator = " =\r\n\t#$,;\"";  //\t is tab key

    std::string name = textFileParser.ReadNextWord( separator );

    int arraySize = 1;
    std::string * value = new std::string[ arraySize ];

    value[ 0 ] = textFileParser.ReadNextWord( separator );

    ONEFLOW::ProcessData( name, value, keyWordIndex, arraySize );

    delete[] value;

    return arraySize;
}

int GetParameterArraySize( const std::string & word )
{
    int arraySize = - 1;
    if ( Word::IsDigit( word ) )
    {
        arraySize = StringToDigit< int >( word );
    }
    else
    {
        arraySize = GetDataValue< int >( word );
    }
    return arraySize;
}

void ReadOneFLOWScriptFile( const std::string & fileName )
{
    TextFileParser textFileParser;

    textFileParser.OpenFile( fileName, std::ios_base::in );

    ONEFLOW::ReadOneFLOWScriptFile( textFileParser );

    textFileParser.CloseFile();
}

void mytestjson();

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

//void ReadOneFLOWScriptFile( const std::string & fileName )
//{
//    std::string jsonFileName = GetJsonFileName( fileName );
//
//    TextFileParser textFileParser;
//    textFileParser.OpenFile( fileName, std::ios_base::in );
//
//    //string name, word;
//
//    //\t is the tab key
//    std::string keyWordSeparator = " =\r\n\t#$,;\"";
//
//    textFileParser.SetDefaultSeparator( keyWordSeparator );

//    DataBaseType::Init();
//
//    Json::Value jsonRoot;
//
//    while ( ! textFileParser.ReachTheEndOfFile() )
//    {
//        bool resultFlag = textFileParser.ReadNextMeaningfulLine();
//        if ( ! resultFlag ) break;
//
//        std::string keyWord = textFileParser.ReadNextWord();
//
//        if ( keyWord == "" ) continue;
//
//        int keyWordIndex = DataBaseType::GetIndex( keyWord );
//
//        Json::Value jsonItem;
//        std::string varName;
//        std::vector< std::string > varArray;
//        GetParaInfo( textFileParser, varName, varArray );
//        jsonItem[ "type" ] = keyWord;
//        jsonItem[ "value" ] = varArray[0];
//        
//        jsonRoot[ varName ] = jsonItem;
//
//        ONEFLOW::ProcessData( varName, &varArray[0], keyWordIndex, varArray.size() );
//    }
//
//    //std::cout << jsonRoot.toStyledString() << std::endl;
//
//    //std::ofstream ofs;
//    //ofs.open( jsonFileName.c_str() );
//    //ofs << jsonRoot.toStyledString();
//    //ofs.close();
//
//    textFileParser.CloseFile();
//}

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
                Stop( errorMessage );
            }
        }
    }
}

void mytestjson()
{
}

void ReadControlInfo()
{
    HXBcastString( Prj::prjBaseDir, Parallel::serverid );

    if ( Parallel::IsServer() )
    {
        ONEFLOW::ReadPrjScript();
    }

    Parallel::TestSayHelloFromEveryProcess();
    ONEFLOW::BroadcastControlParameterToAllProcessors();
    ONEFLOW::DumpDataBase();
}

void DumpDataBase()
{
    DataBase * dataBase = ONEFLOW::GetGlobalDataBase();
    std::fstream file;
    std::string fileName = "/log/database.log";
    PIO::OpenPrjFile( file, fileName, std::ios_base::out );
    dataBase->dataPara->DumpData( file );
    PIO::CloseFile( file );
}

void ReadPrjScript()
{
    std::vector< std::string > scriptFileNameList;
    ONEFLOW::ReadScriptFileNameList( scriptFileNameList );
    ONEFLOW::ReadMultiScriptFiles( scriptFileNameList );
}

void ReadScriptFileNameList( std::vector< std::string > & scriptFileNameList )
{
    TextFileParser textFileParser;

    textFileParser.OpenPrjFile( "script/control.txt", std::ios_base::in );

    //\t is Tab Key
    std::string keyWordSeparator = " ()\r\n\t#$,;\"";
    textFileParser.SetDefaultSeparator( keyWordSeparator );

    while ( ! textFileParser.ReachTheEndOfFile()  )
    {
        bool flag = textFileParser.ReadNextNonEmptyLine();
        if ( ! flag ) break;
        std::string scriptFileName = textFileParser.ReadNextWord();
        OStream &logger = OStream::Instance();
        logger.ClearAll();
        logger << Prj::prjBaseDir << "script/" << scriptFileName;
        std::string fullScriptFileName = logger.str();
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

void CompressData( DataBook *& dataBook )
{
    DataBase * globalDataBase = ONEFLOW::GetGlobalDataBase();

    ONEFLOW::CompressData( globalDataBase, dataBook );
}

void DecompressData( DataBook * dataBook )
{
    DataBase * globalDataBase = ONEFLOW::GetGlobalDataBase();
    ONEFLOW::DecompressData( globalDataBase, dataBook );
}

void CompressData( DataBase * dataBase, DataBook *& dataBook )
{
    DataPara::DataSET * dataSet = dataBase->dataPara->GetDataSet();
    DataPara::DataSET::iterator iter;

    int ndata = static_cast<int> (dataSet->size());

    ONEFLOW::HXWrite( dataBook, ndata );

    for ( iter = dataSet->begin(); iter != dataSet->end(); ++ iter )
    {
        DataV * datav = ( * iter );
        ONEFLOW::HXWriteDataV( dataBook, datav );
    }
}

void DecompressData( DataBase * dataBase, DataBook * dataBook )
{
    DataPara::DataSET * dataSet = dataBase->dataPara->GetDataSet();
    DataPara::DataSET::iterator iter;

    dataBook->MoveToBegin();

    int ndata = 0;
    ONEFLOW::HXRead( dataBook, ndata );

    for ( int i = 0; i < ndata; ++ i )
    {
        DataV * datav = new DataV();
        ONEFLOW::HXReadDataV( dataBook, datav );
        dataBase->dataPara->UpdateDataPointer( datav );
    }
}

EndNameSpace

/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2019 He Xin and the OneFLOW contributors.
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
#include "BasicIO.h"
#include "Prj.h"
#include <iostream>
#include <string>
#include <map>
using namespace std;

BeginNameSpace( ONEFLOW )

bool IsArrayParameter( const string & lineOfName )
{
    const string::size_type npos = - 1;

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

void ReadBasicData( AsciiFileRead & asciiFileRead )
{
    //string name, word;

    //\tÎªtab¼ü
    string keyWordSeparator = " =\r\n\t#$,;\"";

    asciiFileRead.SetDefaultSeparator( keyWordSeparator );

    map < string, int > keyWordMap;

    keyWordMap.insert( pair< string, int >( "int", HX_INT ) );
    keyWordMap.insert( pair< string, int >( "float", HX_FLOAT ) );
    keyWordMap.insert( pair< string, int >( "double", HX_DOUBLE ) );
    keyWordMap.insert( pair< string, int >( "Real", HX_REAL ) );
    keyWordMap.insert( pair< string, int >( "string", HX_STRING ) );
    keyWordMap.insert( pair< string, int >( "bool", HX_BOOL ) );

    while ( ! asciiFileRead.ReachTheEndOfFile() )
    {
        bool resultFlag = asciiFileRead.ReadNextMeaningfulLine();
        if ( ! resultFlag ) break;

        string keyWord = asciiFileRead.ReadNextWord();

        if ( keyWord == "" ) continue;

        int keyWordIndex = keyWordMap[ keyWord ];

        if ( ONEFLOW::IsArrayParameter( asciiFileRead.GetCurrentLine() ) )
        {
            ONEFLOW::AnalysisArrayParameter( asciiFileRead, keyWordIndex );
        }
        else
        {
            ONEFLOW::AnalysisScalarParameter( asciiFileRead, keyWordIndex );
        }
    }
}

void AnalysisArrayParameter( AsciiFileRead & asciiFileRead, int keyWordIndex )
{
    string errorMessage = "error in parameter file";
    string arrayParameterSeparator = "=\r\n\t#$,;\"";

    string compositeArrayParameterNameInformation = asciiFileRead.ReadNextWord( arrayParameterSeparator );

    //Array pattern
    string compositeArrayParameterSeparator = " =\r\n\t#$,;\"[]";
    string arrayParameterName, arraySizeInformation;

    arrayParameterName = ONEFLOW::FindNextWord( compositeArrayParameterNameInformation, compositeArrayParameterSeparator );
    arraySizeInformation = ONEFLOW::FindNextWord( compositeArrayParameterNameInformation, compositeArrayParameterSeparator );

    int arraySize = ONEFLOW::GetParameterArraySize( arraySizeInformation );

    string * valueContainer = new string[ arraySize ];

    for ( int i = 0; i < arraySize; ++ i )
    {
        valueContainer[ i ] = asciiFileRead.ReadNextWord( compositeArrayParameterSeparator );
        //It shows that these contents can't be written within 1 lines
        if ( valueContainer[ i ] == "" )
        {
            asciiFileRead.ReadNextNonEmptyLine();
            valueContainer[ i ] = asciiFileRead.ReadNextWord( compositeArrayParameterSeparator );
            if ( valueContainer[ i ] == "" )
            {
                Stop( errorMessage );
            }
        }
    }
    ONEFLOW::ProcessData( arrayParameterName, valueContainer, keyWordIndex, arraySize );

    delete[] valueContainer;
}

int AnalysisScalarParameter( AsciiFileRead & asciiFileRead, int keyWordIndex )
{
    string errorMessage = "error in parameter file";
    string separator = " =\r\n\t#$,;\"";  //\tÎªtab¼ü

    string name = asciiFileRead.ReadNextWord( separator );

    int arraySize = 1;
    string * value = new string[ arraySize ];

    value[ 0 ] = asciiFileRead.ReadNextWord( separator );

    ONEFLOW::ProcessData( name, value, keyWordIndex, arraySize );

    delete[] value;

    return arraySize;
}

int GetParameterArraySize( const string & word )
{
    int arraySize = - 1;
    if ( ONEFLOW::IsDigit( word ) )
    {
        arraySize = StringToDigit< int >( word );
    }
    else
    {
        arraySize = GetDataValue< int >( word );
    }
    return arraySize;
}

void ReadHXFile( const std::string & fileName )
{
    AsciiFileRead asciiFileRead;

    asciiFileRead.OpenFile( fileName, ios_base::in );

    ONEFLOW::ReadBasicData( asciiFileRead );

    asciiFileRead.CloseFile();
}

void ReadControlInfo()
{
    if ( Parallel::IsServer() )
    {
        ONEFLOW::ReadPrjBaseDir();
    }

    HXBcastString( PrjStatus::prjBaseDir, Parallel::serverid );

    if ( Parallel::IsServer() )
    {
        ONEFLOW::ReadHXScript();
    }

    Parallel::TestSayHelloFromEveryProcess();
    ONEFLOW::BroadcastControlParameterToAllProcessors();
}

void ReadPrjBaseDir()
{
    if ( PrjStatus::prjBaseDir != "" ) return;

    string baseDir = "./";

    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << baseDir << "/prjFile.txt";

    string prjFile = ONEFLOW::StrIO.str();

    ONEFLOW::ReadHXFile( prjFile );

    string prjName = ONEFLOW::GetDataValue< string >( "prjName" );

    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << "./" << baseDir << "/" << prjName << "/";

    PrjStatus::prjBaseDir = ONEFLOW::StrIO.str();
}

void ReadHXScript()
{
    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << PrjStatus::prjBaseDir << "script/control.txt";

    string fileName = ONEFLOW::StrIO.str();
    ONEFLOW::ReadHXFile( fileName );
    ONEFLOW::ReadMultiFile();
}

int GetNumberOfParameterFiles()
{
    return ONEFLOW::GetDataValue< int >( "numberOfParameterFiles" );
}

void ReadMultiFile()
{
    int numberOfParameterFiles = ONEFLOW::GetNumberOfParameterFiles();

    for ( int iFile = 0; iFile < numberOfParameterFiles; ++ iFile )
    {
        string fileName = ONEFLOW::GetParameterFileName( iFile );

        ONEFLOW::StrIO.ClearAll();
        ONEFLOW::StrIO << PrjStatus::prjBaseDir << fileName;

        fileName = ONEFLOW::StrIO.str();

        ONEFLOW::ReadHXFile( fileName );
    }
}

std::string GetParameterFileName( int iFile )
{
    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << "parameterFileName" << iFile;
    std::string fileNameString = ONEFLOW::StrIO.str();

    return ONEFLOW::GetDataValue< std::string >( fileNameString );
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
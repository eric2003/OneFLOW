/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
#include <iostream>
#include <string>
#include <vector>
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

void ReadOneFLOWScriptFile( FileIO & fileIO )
{
    //string name, word;

    //\tÎªtab¼ü
    string keyWordSeparator = " =\r\n\t#$,;\"";

    fileIO.SetDefaultSeparator( keyWordSeparator );

    map < string, int > keyWordMap;

    keyWordMap.insert( pair< string, int >( "int", HX_INT ) );
    keyWordMap.insert( pair< string, int >( "float", HX_FLOAT ) );
    keyWordMap.insert( pair< string, int >( "double", HX_DOUBLE ) );
    keyWordMap.insert( pair< string, int >( "Real", HX_REAL ) );
    keyWordMap.insert( pair< string, int >( "string", HX_STRING ) );
    keyWordMap.insert( pair< string, int >( "bool", HX_BOOL ) );

    while ( ! fileIO.ReachTheEndOfFile() )
    {
        bool resultFlag = fileIO.ReadNextMeaningfulLine();
        if ( ! resultFlag ) break;

        string keyWord = fileIO.ReadNextWord();

        if ( keyWord == "" ) continue;

        int keyWordIndex = keyWordMap[ keyWord ];

        if ( ONEFLOW::IsArrayParameter( fileIO.GetCurrentLine() ) )
        {
            ONEFLOW::AnalysisArrayParameter( fileIO, keyWordIndex );
        }
        else
        {
            ONEFLOW::AnalysisScalarParameter( fileIO, keyWordIndex );
        }
    }
}

void AnalysisArrayParameter( FileIO & fileIO, int keyWordIndex )
{
    string errorMessage = "error in parameter file";
    string commSeparator = "=\r\n\t#$,;\"";

    string ayrrayInfo = fileIO.ReadNextWord( commSeparator );

    //Array pattern
    string arraySeparator = " =\r\n\t#$,;\"[]";
    string arrayName, arraySizeName;

    arrayName = Word::FindNextWord( ayrrayInfo, arraySeparator );
    arraySizeName = Word::FindNextWord( ayrrayInfo, arraySeparator );

    int arraySize = ONEFLOW::GetParameterArraySize( arraySizeName );

    string * valueContainer = new string[ arraySize ];

    for ( int i = 0; i < arraySize; ++ i )
    {
        valueContainer[ i ] = fileIO.ReadNextWord( arraySeparator );
        //It shows that these contents can't be written within 1 lines
        if ( valueContainer[ i ] == "" )
        {
            fileIO.ReadNextNonEmptyLine();
            valueContainer[ i ] = fileIO.ReadNextWord( arraySeparator );
            if ( valueContainer[ i ] == "" )
            {
                Stop( errorMessage );
            }
        }
    }
    ONEFLOW::ProcessData( arrayName, valueContainer, keyWordIndex, arraySize );

    delete[] valueContainer;
}

int AnalysisScalarParameter( FileIO & fileIO, int keyWordIndex )
{
    string errorMessage = "error in parameter file";
    string separator = " =\r\n\t#$,;\"";  //\t is tab key

    string name = fileIO.ReadNextWord( separator );

    int arraySize = 1;
    string * value = new string[ arraySize ];

    value[ 0 ] = fileIO.ReadNextWord( separator );

    ONEFLOW::ProcessData( name, value, keyWordIndex, arraySize );

    delete[] value;

    return arraySize;
}

int GetParameterArraySize( const string & word )
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
    FileIO fileIO;

    fileIO.OpenFile( fileName, ios_base::in );

    ONEFLOW::ReadOneFLOWScriptFile( fileIO );

    fileIO.CloseFile();
}

void ReadControlInfo()
{
    HXBcastString( PrjStatus::prjBaseDir, Parallel::serverid );

    if ( Parallel::IsServer() )
    {
        ONEFLOW::ReadPrjScript();
    }

    Parallel::TestSayHelloFromEveryProcess();
    ONEFLOW::BroadcastControlParameterToAllProcessors();
}

void ReadPrjScript()
{
    vector< string > scriptFileNameList;
    ONEFLOW::ReadScriptFileNameList( scriptFileNameList );
    ONEFLOW::ReadMultiScriptFiles( scriptFileNameList );
}

void ReadScriptFileNameList( vector< string > & scriptFileNameList )
{
    FileIO ioFile;

    ioFile.OpenPrjFile( "script/control.txt", ios_base::in );

    //\t is Tab Key
    string keyWordSeparator = " ()\r\n\t#$,;\"";
    ioFile.SetDefaultSeparator( keyWordSeparator );

    while ( ! ioFile.ReachTheEndOfFile()  )
    {
        bool flag = ioFile.ReadNextNonEmptyLine();
        if ( ! flag ) break;
        string scriptFileName = ioFile.ReadNextWord();
        ONEFLOW::StrIO.ClearAll();
        ONEFLOW::StrIO << PrjStatus::prjBaseDir << "script/" << scriptFileName;
        string fullScriptFileName = ONEFLOW::StrIO.str();
        scriptFileNameList.push_back( fullScriptFileName );
    }

    ioFile.CloseFile();
}

void ReadMultiScriptFiles( vector< string > & scriptFileNameList )
{
    int numberOfParameterFiles = scriptFileNameList.size();

    for ( int iFile = 0; iFile < numberOfParameterFiles; ++ iFile )
    {
        string & scriptFileName = scriptFileNameList[ iFile ];

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
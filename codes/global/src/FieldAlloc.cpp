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
#include "FieldAlloc.h"
#include "Prj.h"
#include "Fatal.h"
#include "FieldImp.h"
#include "UsdPara.h"
#include "SolverInfo.h"
#include "SolverDef.h"
#include "TextFileParser.h"
#include "OStream.h"
#include "DataBase.h"
#include "RegisterUtils.h"
#include "Zone.h"
#include "Grid.h"
#include "InterFace.h"

BeginNameSpace( ONEFLOW )

void FieldAlloc::AllocateAllFields( int solverType, const std::string & basicString )
{
    FieldAlloc::RegisterInterfaceVar( solverType, basicString );
    FieldAlloc::AllocateGlobalField( solverType, basicString );
    FieldAlloc::InitField( solverType, basicString );
}

void FieldAlloc::InitField( int solverType, const std::string & basicString )
{
    std::string fileName = Prj::GetSystemFileName( basicString + "/alloc/init.txt" );
    BoolIO boolIO;
    boolIO.ReadFile( fileName, 1 );

    SetFieldValues(
        solverType,
        boolIO.nameValuePair );
}


void FieldAlloc::RegisterInterfaceVar( int solverType, const std::string & basicString )
{
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( solverType );
    if ( solverInfo->registerInterface ) return;
    solverInfo->registerInterface = 1;

    StringField fileNameList;
    IntField fieldTypeList;

    FieldAlloc::CalcInterfaceFileName( basicString, fileNameList );
    FieldAlloc::CalcInterfaceFileType( fieldTypeList );
    
    for ( int iFile = 0; iFile < fileNameList.size(); ++ iFile )
    {
        BoolIO boolIO;
        boolIO.ReadFile( fileNameList[ iFile ] );
        int fieldType = fieldTypeList[ iFile ];

        AddInterfaceFieldNames(
            solverType,
            fieldType,
            boolIO.nameValuePair.nameList );
    }
}

void FieldAlloc::AllocateGlobalField( int solverType, const std::string & basicString )
{
    FieldFactory::AddFieldManager( solverType );
    StringField fileNameList;
    FieldAlloc::CalcInnerFieldFileName( basicString, fileNameList );

    for ( int iFile = 0; iFile < fileNameList.size(); ++ iFile )
    {
        ReadSuperPara readSuperPara;

        readSuperPara.solverType = solverType;
        readSuperPara.Register( fileNameList[ iFile ], iFile );
    }

    FieldAlloc::AllocateAllKindsOfInterfaceField( solverType );
}

void FieldAlloc::AllocateAllKindsOfInterfaceField( int solverType )
{
    FieldManager * fieldManager = FieldFactory::GetFieldManager( solverType );
    fieldManager->AllocateInnerAndBcField();
    FieldAlloc::AllocateInterfaceField( fieldManager->iFieldProperty.get() );
    FieldAlloc::AllocateOversetInterfaceField( fieldManager->iFieldProperty.get() );
}

void FieldAlloc::AllocateInterfaceField( IFieldProperty * iFieldProperty )
{
    Grid * grid = Zone::GetGrid();

    InterFace * interFace = grid->interFace;

    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    int nIFaces = grid->interFace->nIFaces;
    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        iFieldProperty->AllocateInterfaceField( nIFaces, interFace->dataSend[ ghostId ] );
        iFieldProperty->AllocateInterfaceField( nIFaces, interFace->dataRecv[ ghostId ] );
    }
}

void FieldAlloc::AllocateOversetInterfaceField( IFieldProperty * iFieldProperty )
{
}

void FieldAlloc::CalcInnerFieldFileName( const std::string & basicString, StringField & fileNameList )
{
    StringField basicNameList;
    basicNameList.push_back( "unsteady" );
    basicNameList.push_back( "inner"    );
    basicNameList.push_back( "face"     );
    basicNameList.push_back( "bc"       );

    //std::string rootString = logger.str();

    std::string rootString = Prj::GetSystemFileName( basicString + "/alloc/" );

    OStream &logger = OStream::Instance();

    for ( int i = 0; i < basicNameList.size(); ++ i )
    {
        logger.ClearAll();
        logger << rootString << basicNameList[ i ] << ".txt";

        std::string name = logger.str();

        fileNameList.push_back( name );
    }
}

void FieldAlloc::CalcInterfaceFileName( const std::string & basicString, StringField & fileNameList )
{
    StringField basicNameList;
    basicNameList.push_back( "inter"        );
    basicNameList.push_back( "interDq"      );
    basicNameList.push_back( "interGrad"    );
    basicNameList.push_back( "interOverset" );

    std::string rootString = Prj::GetSystemFileName( basicString + "/alloc/" );
    OStream &logger = OStream::Instance();

    for ( int i = 0; i < basicNameList.size(); ++ i )
    {
        logger.ClearAll();
        logger << rootString << basicNameList[ i ] << ".txt";

        std::string name = logger.str();

        fileNameList.push_back( name );
    }
}

void FieldAlloc::CalcInterfaceFileType( IntField & fieldTypeList )
{
    fieldTypeList.push_back( ONEFLOW::INTERFACE_DATA          );
    fieldTypeList.push_back( ONEFLOW::INTERFACE_DQ_DATA       );
    fieldTypeList.push_back( ONEFLOW::INTERFACE_GRADIENT_DATA );
    fieldTypeList.push_back( ONEFLOW::INTERFACE_OVERSET_DATA  );
}


void SetFieldValues(
    int solverType,
    const NameValuePair & valuePair )
{
    FieldManager * fieldManager =
        FieldFactory::GetFieldManager( solverType );

    int nVar = valuePair.nameList.size();

    for ( int iVar = 0; iVar < nVar; ++ iVar )
    {
        fieldManager->SetField(
            valuePair.nameList[ iVar ],
            valuePair.valueList[ iVar ] );
    }
}


bool CalcBoolExp( bool var1, const std::string & opName, bool var2 )
{
    if ( opName == "&&" )
    {
        return var1 && var2;
    }
    else if ( opName == "||" )
    {
        return var1 || var2;
    }
    return false;
}

bool CalcBoolExp( const std::string & varName1, const std::string & opName, const std::string & varName2 )
{
    int var1 = ONEFLOW::GetVarDimension( varName1 );
    int var2 = ONEFLOW::GetVarDimension( varName2 );
    if ( opName == ">" )
    {
        return var1 > var2;
    }
    else if ( opName == ">=" )
    {
        return var1 >= var2;
    }
    else if ( opName == "==" )
    {
        return var1 == var2;
    }
    else if ( opName == "<" )
    {
        return var1 < var2;
    }
    else if ( opName == "<=" )
    {
        return var1 <= var2;
    }
    else if ( opName == "!=" )
    {
        return var1 != var2;
    }
    return false;
}

bool CalcVarValue( const std::string & varName, StringField & boolName, BoolField & boolVar )
{
    int index = -1;
    for ( int i = 0; i < boolName.size(); ++ i )
    {
        if ( varName == boolName[ i ] )
        {
            index = i;
            break;
        }
    }
    return boolVar[ index ];
}

int GetVarDimension( const std::string & dimName )
{
    if ( Word::IsDigit( dimName ) )
    {
        return StringToDigit< int >( dimName );
    }
    else
    {
        return GetDataValue< int >( dimName );
    }
}

BoolIO::BoolIO()
{
    ;
}

BoolIO::~BoolIO()
{
    ;
}

void BoolIO::Add( const std::string & name, bool value )
{
    this->boolNameList.push_back( name );
    this->boolValueList.push_back( value );
}

void BoolIO::ReadBool( TextFileParser & textFileParser )
{
    std::string varName = textFileParser.ReadNextWord();
    std::string word    = textFileParser.ReadNextWord();
    std::string var1    = textFileParser.ReadNextWord();
    std::string opName  = textFileParser.ReadNextWord();
    std::string var2    = textFileParser.ReadNextWord();

    bool boolValue = ONEFLOW::CalcBoolExp( var1, opName, var2 );

    this->Add( varName, boolValue );
}


void BoolIO::ReadSuperBool( TextFileParser & textFileParser )
{
    std::string varName = textFileParser.ReadNextWord();
    std::string word    = textFileParser.ReadNextWord();
    std::string var1    = textFileParser.ReadNextWord();
    std::string opName  = textFileParser.ReadNextWord();
    std::string var2    = textFileParser.ReadNextWord();

    bool varVaule1 =
        ONEFLOW::CalcVarValue(
            var1,
            this->boolNameList,
            this->boolValueList );

    bool varVaule2 =
        ONEFLOW::CalcVarValue(
            var2,
            this->boolNameList,
            this->boolValueList );

    bool boolValue =
        ONEFLOW::CalcBoolExp(
            varVaule1,
            opName,
            varVaule2 );

    this->Add( varName, boolValue );
}

bool BoolIO::CalcVarValue( const std::string & varName )
{
    bool result = ONEFLOW::CalcVarValue( varName, this->boolNameList, this->boolValueList );
    return result;
}

void BoolIO::Read(
    TextFileParser & textFileParser,
    int valueFlag,
    ParaNameDimData * paraNameDimData )
{
    if ( valueFlag == 0 )
    {
        std::string varName = textFileParser.ReadNextWord();
        nameValuePair.nameList.push_back( varName );
    }
    else if ( valueFlag == 1 )
    {
        std::string varName = textFileParser.ReadNextWord();
        nameValuePair.nameList.push_back( varName );

        Real varValue = textFileParser.ReadNextDigit< Real >();
        nameValuePair.valueList.push_back( varValue );
    }
    else if ( valueFlag == 2 )
    {
        if ( paraNameDimData == nullptr )
        {
            Fatal( "ParaNameDimData is required when BoolIO valueFlag is 2." );
        }

        std::string varName      = textFileParser.ReadNextWord();
        std::string varDimension = textFileParser.ReadNextWord();
        std::string typeName     = textFileParser.ReadNextWord();

        int dimension = ONEFLOW::GetVarDimension( varDimension );

        ParaNameDim * paraNameDim =
            paraNameDimData->GetParaNameDim( typeName );

        paraNameDim->nameList.push_back( varName );
        paraNameDim->dimList.push_back( dimension );
    }
}

void BoolIO::ReadFile(
    const std::string & fileName,
    int valueFlag,
    ParaNameDimData * paraNameDimData )
{
    // \t is the tab key
    std::string separator = " \r\n\t#$,;\"()";

    TextFileParser textFileParser;
    textFileParser.OpenFile( fileName, std::ios_base::in );
    textFileParser.SetDefaultSeparator( separator );

    while ( ! textFileParser.ReachTheEndOfFile() )
    {
        bool flag = textFileParser.ReadNextNonEmptyLine();
        if ( ! flag ) break;

        std::string keyWord = textFileParser.ReadNextWord();

        if ( keyWord == "true" )
        {
            this->Read(
                textFileParser,
                valueFlag,
                paraNameDimData );
        }
        else if ( keyWord == "bool" )
        {
            this->ReadBool( textFileParser );
        }
        else if ( keyWord == "superbool" )
        {
            this->ReadSuperBool( textFileParser );
        }
        else
        {
            std::string expression = keyWord;
            bool flag = this->CalcVarValue( expression );

            if ( flag )
            {
                this->Read(
                    textFileParser,
                    valueFlag,
                    paraNameDimData );
            }
        }
    }

    textFileParser.CloseFile();
}

void AddInterfaceFieldNames(
    int solverType,
    int fieldType,
    const StringField & nameList )
{
    VarNameSolver * varNameSolver =
        VarNameFactory::GetVarNameSolver(
            solverType,
            fieldType );

    int numberOfVariables = nameList.size();

    for ( int iVariable = 0;
        iVariable < numberOfVariables;
        ++ iVariable )
    {
        const std::string & varName =
            nameList[ iVariable ];

        varNameSolver->AddFieldName( varName );
    }
}


ParaNameDim * ParaNameDimData::GetParaNameDim(
    const std::string & typeName )
{
    if ( typeName == "all" )
    {
        return &this->comPara;
    }
    else if ( typeName == "str" )
    {
        return &this->strPara;
    }
    else
    {
        return &this->unsPara;
    }
}

ReadSuperPara::ReadSuperPara()
    : paraNameDimData( std::make_unique<ParaNameDimData>() )
{
}

ReadSuperPara::~ReadSuperPara() = default;

void ReadSuperPara::AddUnsteadyInnerFieldProperty()
{
    this->AddFieldProperties( FieldLocation::Inner );
    FieldManager * fieldManager = FieldFactory::GetFieldManager( this->solverType );

    UsdPara * usdPara = fieldManager->usdPara.get();
    int nEqu = this->paraNameDimData->comPara.dimList[ 0 ];
    usdPara->Init( this->paraNameDimData->comPara.nameList, nEqu );
}

void ReadSuperPara::AddFieldProperties(
    FieldLocation location )
{
    this->AddBasicFieldProperty(
        &this->paraNameDimData->unsPara,
        location,
        FieldCategory::Unstructured );

    this->AddBasicFieldProperty(
        &this->paraNameDimData->strPara,
        location,
        FieldCategory::Structured );

    this->AddBasicFieldProperty(
        &this->paraNameDimData->comPara,
        location,
        FieldCategory::Common );
}

void ReadSuperPara::AddBasicFieldProperty(
    ParaNameDim * paraNameDim,
    FieldLocation location,
    FieldCategory category )
{
    FieldManager * fieldManager =
        FieldFactory::GetFieldManager( this->solverType );

    int nVar = paraNameDim->nameList.size();

    for ( int iVar = 0; iVar < nVar; ++ iVar )
    {
        const std::string & varName =
            paraNameDim->nameList[ iVar ];

        int nEqu =
            paraNameDim->dimList[ iVar ];

        fieldManager->AddField(
            varName,
            nEqu,
            category,
            location );
    }
}

void ReadSuperPara::Register( const std::string & fileName, int index )
{
    BoolIO boolIO;

    boolIO.ReadFile(
        fileName,
        2,
        this->paraNameDimData.get() );

    switch ( index )
    {
    case 0:
        // Unsteady fields have separate semantics.
        this->AddUnsteadyInnerFieldProperty();
        break;

    case 1:
        this->AddFieldProperties( FieldLocation::Inner );
        break;

    case 2:
        this->AddFieldProperties( FieldLocation::Face );
        break;

    case 3:
        this->AddFieldProperties( FieldLocation::Boundary );
        break;
    }
}


EndNameSpace

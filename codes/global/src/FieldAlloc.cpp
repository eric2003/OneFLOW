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
#include "FieldAlloc.h"
#include "SimuCtrl.h"
#include "FieldImp.h"
#include "UsdPara.h"
#include "SolverInfo.h"
#include "SolverDef.h"
#include "FileIO.h"
#include "OStream.h"
#include "DataBase.h"
#include "RegisterUtil.h"
#include "Zone.h"
#include "Grid.h"
#include "InterFace.h"

BeginNameSpace( ONEFLOW )

FieldAlloc::FieldAlloc()
{
    ;
}

FieldAlloc::~FieldAlloc()
{
    ;
}

void FieldAlloc::AllocateAllFields( int sTid, const string & basicString )
{
    FieldAlloc::RegisterInterfaceVar( sTid, basicString );
    FieldAlloc::AllocateGlobalField( sTid, basicString );
    FieldAlloc::InitField( sTid, basicString );
}

void FieldAlloc::InitField( int sTid, const string & basicString )
{
    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << SimuCtrl::system_root << basicString << "/alloc/" << "init.txt";
    std::string fileName = ONEFLOW::StrIO.str();

    BoolIO boolIO;
    boolIO.ReadFile( fileName, 1 );

    //FieldNamePair::Read( fileName, nameValuePair );
    FieldNamePair::SetField( sTid, boolIO.nameValuePair );
}


void FieldAlloc::RegisterInterfaceVar( int sTid, const string & basicString )
{
    SolverInfo * solverInfo = SolverInfoFactory::GetSolverInfo( sTid );
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

        ReadInterfaceVar::AddFieldName( sTid, fieldType, boolIO.nameValuePair.nameList );
    }
}

void FieldAlloc::AllocateGlobalField( int sTid, const string & basicString )
{
    FieldFactory::AddFieldManager( sTid );
    StringField fileNameList;
    FieldAlloc::CalcInnerFieldFileName( basicString, fileNameList );

    for ( int iFile = 0; iFile < fileNameList.size(); ++ iFile )
    {
        ReadSuperPara * readSuperPara = new ReadSuperPara();

        readSuperPara->sTid = sTid;
        readSuperPara->Register( fileNameList[ iFile ], iFile );
        delete readSuperPara;
    }

    FieldAlloc::AllocateAllKindsOfInterfaceField( sTid );
}

void FieldAlloc::AllocateAllKindsOfInterfaceField( int sTid )
{
    FieldManager * fieldManager = FieldFactory::GetFieldManager( sTid );
    fieldManager->AllocateInnerAndBcField();
    FieldAlloc::AllocateInterfaceField( fieldManager->iFieldProperty );
    FieldAlloc::AllocateOversetInterfaceField( fieldManager->iFieldProperty );
}

void FieldAlloc::AllocateInterfaceField( IFieldProperty * iFieldProperty )
{
    Grid * grid = Zone::GetGrid();

    InterFace * interFace = grid->interFace;

    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    int nIFace = grid->interFace->nIFace;
    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        iFieldProperty->AllocateInterfaceField( nIFace, interFace->dataSend[ ghostId ] );
        iFieldProperty->AllocateInterfaceField( nIFace, interFace->dataRecv[ ghostId ] );
    }
}

void FieldAlloc::AllocateOversetInterfaceField( IFieldProperty * iFieldProperty )
{
}

void FieldAlloc::CalcInnerFieldFileName( const string & basicString, StringField & fileNameList )
{
    StringField basicNameList;
    basicNameList.push_back( "unsteady" );
    basicNameList.push_back( "inner"    );
    basicNameList.push_back( "face"     );
    basicNameList.push_back( "bc"       );

    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << SimuCtrl::system_root << basicString << "/alloc/";
    string rootString = ONEFLOW::StrIO.str();

    for ( int i = 0; i < basicNameList.size(); ++ i )
    {
        ONEFLOW::StrIO.ClearAll();
        ONEFLOW::StrIO << rootString << basicNameList[ i ] << ".txt";

        string name = ONEFLOW::StrIO.str();

        fileNameList.push_back( name );
    }
}

void FieldAlloc::CalcInterfaceFileName( const string & basicString, StringField & fileNameList )
{
    StringField basicNameList;
    basicNameList.push_back( "inter"        );
    basicNameList.push_back( "interDq"      );
    basicNameList.push_back( "interGrad"    );
    basicNameList.push_back( "interOverset" );

    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << SimuCtrl::system_root << basicString << "/alloc/";
    string rootString = ONEFLOW::StrIO.str();

    for ( int i = 0; i < basicNameList.size(); ++ i )
    {
        ONEFLOW::StrIO.ClearAll();
        ONEFLOW::StrIO << rootString << basicNameList[ i ] << ".txt";

        string name = ONEFLOW::StrIO.str();

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

FieldNamePair::FieldNamePair()
{
}

FieldNamePair::~FieldNamePair()
{
}

void FieldNamePair::SetField( int sTid, NameValuePair & valuePair )
{
    FieldManager * fieldManager = FieldFactory::GetFieldManager( sTid );
    int nVar = valuePair.nameList.size();
    for ( int iVar = 0; iVar < nVar; ++ iVar )
    {
        fieldManager->SetField( valuePair.nameList[ iVar ], valuePair.valueList[ iVar ] );
    }
}

bool CalcBoolExp( bool var1, const string & opName, bool var2 )
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

bool CalcBoolExp( const string & varName1, const string & opName, const string & varName2 )
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

bool CalcVarValue( const string & varName, StringField & boolName, BoolField & boolVar )
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

int GetVarDimension( const string & dimName )
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
    this->valueFlag = 0;
}

BoolIO::~BoolIO()
{
}

void BoolIO::Add( const string & name, bool value )
{
    this->boolNameList.push_back( name );
    this->boolValueList.push_back( value );
}

void BoolIO::ReadBool( FileIO * ioFile )
{
    string varName = ioFile->ReadNextWord();
    string word    = ioFile->ReadNextWord(); //"="
    string var1    = ioFile->ReadNextWord();
    string opName  = ioFile->ReadNextWord();
    string var2    = ioFile->ReadNextWord();

    bool boolValue = ONEFLOW::CalcBoolExp( var1, opName, var2 );

    this->Add( varName, boolValue );
}

void BoolIO::ReadSuperBool( FileIO * ioFile )
{
    string varName = ioFile->ReadNextWord();
    string word    = ioFile->ReadNextWord(); //"="
    string var1    = ioFile->ReadNextWord();
    string opName  = ioFile->ReadNextWord();
    string var2    = ioFile->ReadNextWord();

    bool varVaule1 = ONEFLOW::CalcVarValue( var1, this->boolNameList, this->boolValueList );
    bool varVaule2 = ONEFLOW::CalcVarValue( var2, this->boolNameList, this->boolValueList );

    bool boolValue = ONEFLOW::CalcBoolExp( varVaule1, opName, varVaule2 );

    this->Add( varName, boolValue );
}

bool BoolIO::CalcVarValue( const string & varName )
{
    bool result = ONEFLOW::CalcVarValue( varName, this->boolNameList, this->boolValueList );
    return result;
}

void BoolIO::Read()
{
    if ( valueFlag == 0 )
    {
        string varName = ioFile->ReadNextWord();
        nameValuePair.nameList.push_back( varName );
    }
    else if ( valueFlag == 1 )
    {
        string varName = ioFile->ReadNextWord();
        nameValuePair.nameList.push_back( varName );

        Real varValue = ioFile->ReadNextDigit< Real >();
        nameValuePair.valueList.push_back( varValue );
    }
    else if ( valueFlag == 2 )
    {
        string varName      = ioFile->ReadNextWord();
        string varDimension = ioFile->ReadNextWord();
        string typeName     = ioFile->ReadNextWord();

        int dimension = ONEFLOW::GetVarDimension( varDimension );

        ParaNameDim * paraNameDim = paraNameDimData->GetParaNameDim( typeName );
        paraNameDim->nameList.push_back( varName );
        paraNameDim->dimList.push_back( dimension );
    }

}

void BoolIO::ReadFile( const string & fileName, int valueFlag )
{
   //\tÎªtab¼ü
    //string separator  = " =\r\n\t#$,;\"()";
    string separator  = " \r\n\t#$,;\"()";

    FileIO ioFile;
    ioFile.OpenFile( fileName, ios_base::in );
    ioFile.SetDefaultSeparator( separator );

    this->ioFile = & ioFile;
    this->valueFlag = valueFlag;
    while ( ! ioFile.ReachTheEndOfFile()  )
    {
        bool flag = ioFile.ReadNextNonEmptyLine();
        if ( ! flag ) break;
        string keyWord = ioFile.ReadNextWord();

        if ( keyWord == "true" )
        {
            this->Read();
        }
        else if ( keyWord == "bool" )
        {
            this->ReadBool( & ioFile );
        }
        else if ( keyWord == "superbool" )
        {
            this->ReadSuperBool( & ioFile );
        }
        else
        {
            string expression = keyWord;
            bool flag = this->CalcVarValue( expression );

            if ( flag )
            {
                this->Read();
            }
        }
    }

    ioFile.CloseFile();
}

ReadInterfaceVar::ReadInterfaceVar()
{
    ;
}

ReadInterfaceVar::~ReadInterfaceVar()
{
    ;
}

void ReadInterfaceVar::AddFieldName( int sTid, int fieldType, StringField & nameList )
{
    VarNameSolver * varNameSolver = VarNameFactory::GetVarNameSolver( sTid, fieldType );
    int numberOfVariables = nameList.size();
    for ( int iVariable = 0; iVariable < numberOfVariables; ++ iVariable )
    {
        string & varName = nameList[ iVariable ];
        varNameSolver->AddFieldName( varName );
    }
}

ParaNameDim::ParaNameDim()
{
    ;
}

ParaNameDim::~ParaNameDim()
{
    ;
}

ParaNameDimData::ParaNameDimData()
{
    this->comPara = new ParaNameDim();
    this->strPara = new ParaNameDim();
    this->unsPara = new ParaNameDim();
}

ParaNameDimData::~ParaNameDimData()
{
    delete this->comPara;
    delete this->strPara;
    delete this->unsPara;
}

ParaNameDim * ParaNameDimData::GetParaNameDim( const string & typeName )
{
    if ( typeName == "all" )
    {
        return this->comPara;
    }
    else if ( typeName == "str" )
    {
        return this->strPara;
    }
    else
    {
        return this->unsPara;
    }
}

ReadSuperPara::ReadSuperPara()
{
    this->paraNameDimData = new ParaNameDimData();
}

ReadSuperPara::~ReadSuperPara()
{
    delete this->paraNameDimData;
}


void ReadSuperPara::AddInnerFieldProperty()
{
    this->AddBasicFieldProperty( this->paraNameDimData->unsPara, 0, 0 );
    this->AddBasicFieldProperty( this->paraNameDimData->strPara, 0, 1 );
    this->AddBasicFieldProperty( this->paraNameDimData->comPara, 0, 2 );
}

void ReadSuperPara::AddUnsteadyInnerFieldProperty()
{
    this->AddInnerFieldProperty();
    FieldManager * fieldManager = FieldFactory::GetFieldManager( this->sTid );

    UsdPara * usdPara = fieldManager->usdPara;
    int nEqu = this->paraNameDimData->comPara->dimList[ 0 ];
    usdPara->Init( this->paraNameDimData->comPara->nameList, nEqu );
}

void ReadSuperPara::AddFaceFieldProperty()
{
    this->AddBasicFieldProperty( this->paraNameDimData->unsPara, 1, 0 );
    this->AddBasicFieldProperty( this->paraNameDimData->strPara, 1, 1 );
    this->AddBasicFieldProperty( this->paraNameDimData->comPara, 1, 2 );
}

void ReadSuperPara::AddBoundaryFieldProperty()
{
    this->AddBasicFieldProperty( this->paraNameDimData->unsPara, 2, 0 );
    this->AddBasicFieldProperty( this->paraNameDimData->strPara, 2, 1 );
    this->AddBasicFieldProperty( this->paraNameDimData->comPara, 2, 2 );
}

void ReadSuperPara::AddBasicFieldProperty( ParaNameDim * paraNameDim, int fieldType, int type )
{
    FieldManager * fieldManager = FieldFactory::GetFieldManager( this->sTid );

    int nVar = paraNameDim->nameList.size();
    for ( int iVar = 0; iVar < nVar; ++ iVar )
    {
        string & varName = paraNameDim->nameList[ iVar ];
        int nEqu = paraNameDim->dimList[ iVar ];
        if ( fieldType == 0 )
        {
            fieldManager->AddInnerField( varName, nEqu, type );
        }
        else if ( fieldType == 1 )
        {
            fieldManager->AddFaceField( varName, nEqu, type );
        }
        else
        {
            fieldManager->AddBcField( varName, nEqu, type );
        }
    }
}

void ReadSuperPara::Register( const string & fileName, int index )
{
    BoolIO boolIO;
    boolIO.paraNameDimData = this->paraNameDimData;
    boolIO.ReadFile( fileName, 2 );

    if ( index == 0 )
    {
        this->AddUnsteadyInnerFieldProperty();
    }
    else if ( index == 1 )
    {
        this->AddInnerFieldProperty();
    }
    else if ( index == 2 )
    {
        this->AddFaceFieldProperty();
    }
    else if ( index == 3 )
    {
        this->AddBoundaryFieldProperty();
    }
}


EndNameSpace
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

namespace
{
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

    struct FieldFileSpec
    {
        const char * name;
        FieldLocation location;
        bool isUnsteady;
    };

    struct InterfaceFileSpec
    {
        const char * name;
        int fieldType;
    };

    FieldCategory ParseFieldCategory(
        const std::string & typeName )
    {
        if ( typeName == "all" )
        {
            return FieldCategory::Common;
        }

        if ( typeName == "str" )
        {
            return FieldCategory::Structured;
        }

        return FieldCategory::Unstructured;
    }

    void ReadFieldDefinition(
        TextFileParser & textFileParser,
        ParaNameDimData & paraNameDimData )
    {
        std::string varName =
            textFileParser.ReadNextWord();

        std::string varDimension =
            textFileParser.ReadNextWord();

        std::string typeName =
            textFileParser.ReadNextWord();

        int nEqu =
            GetVarDimension( varDimension );

        FieldCategory category =
            ParseFieldCategory( typeName );

        ParaNameDim * paraNameDim =
            paraNameDimData.GetParaNameDim( category );

        paraNameDim->nameList.push_back( varName );
        paraNameDim->nEquList.push_back( nEqu );
    }

    void AddBasicFieldProperty(
        FieldManager * fieldManager,
        ParaNameDim * paraNameDim,
        FieldLocation location,
        FieldCategory category )
    {
        int nVar = paraNameDim->nameList.size();

        for ( int iVar = 0; iVar < nVar; ++ iVar )
        {
            const std::string & varName =
                paraNameDim->nameList[ iVar ];

            int nEqu =
                paraNameDim->nEquList[ iVar ];

            fieldManager->AddField(
                varName,
                nEqu,
                category,
                location );
        }
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

    bool CalcBoolLogic(
        bool var1,
        const std::string & opName,
        bool var2 )
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
    
    bool CompareVar(
        const std::string & varName1,
        const std::string & opName,
        const std::string & varName2 )
    {
        int var1 =
            GetVarDimension( varName1 );

        int var2 =
            GetVarDimension( varName2 );

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

    void ReadFieldDefinitions(
        const std::string & fileName,
        ParaNameDimData & paraNameDimData )
    {
        TextFileParser textFileParser;

        // \t is the tab key
        std::string separator = " \r\n\t#$,;\"()";

        textFileParser.OpenFile(
            fileName,
            std::ios_base::in );

        textFileParser.SetDefaultSeparator(
            separator );

        while ( ! textFileParser.ReachTheEndOfFile() )
        {
            bool flag = textFileParser.ReadNextNonEmptyLine();
            if ( ! flag ) break;

            std::string keyWord =
                textFileParser.ReadNextWord();

            if ( keyWord == "true" )
            {
                ReadFieldDefinition(
                    textFileParser,
                    paraNameDimData );
            }
        }

        textFileParser.CloseFile();
    }

    using BoolLineReader =
        void ( BoolIO::* )( TextFileParser & );

    void ReadBoolFile(
        BoolIO & boolIO,
        const std::string & fileName,
        BoolLineReader trueReader )
    {
        // \t is the tab key
        std::string separator = " \r\n\t#$,;\"()";

        TextFileParser textFileParser;

        textFileParser.OpenFile(
            fileName,
            std::ios_base::in );

        textFileParser.SetDefaultSeparator(
            separator );

        while ( ! textFileParser.ReachTheEndOfFile() )
        {
            bool flag =
                textFileParser.ReadNextNonEmptyLine();

            if ( ! flag ) break;

            std::string keyWord =
                textFileParser.ReadNextWord();

            if ( keyWord == "true" )
            {
                ( boolIO.*trueReader )(
                    textFileParser );
            }
            else if ( keyWord == "bool" )
            {
                boolIO.ReadBool(
                    textFileParser );
            }
            else if ( keyWord == "superbool" )
            {
                boolIO.ReadSuperBool(
                    textFileParser );
            }
            else
            {
                bool flag =
                    boolIO.GetBoolValue( keyWord );

                if ( flag )
                {
                    ( boolIO.*trueReader )(
                        textFileParser );
                }
            }
        }

        textFileParser.CloseFile();
    }
}

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
    boolIO.ReadValueFile( fileName );

    SetFieldValues(
        solverType,
        boolIO.GetNameValuePair() );
}

void FieldAlloc::RegisterInterfaceVar(
    int solverType,
    const std::string & basicString )
{
    SolverInfo * solverInfo =
        SolverInfoFactory::GetSolverInfo( solverType );

    if ( solverInfo->registerInterface ) return;

    solverInfo->registerInterface = 1;

    const InterfaceFileSpec interfaceFileSpecs[] =
    {
        { "inter",        ONEFLOW::INTERFACE_DATA          },
        { "interDq",      ONEFLOW::INTERFACE_DQ_DATA       },
        { "interGrad",    ONEFLOW::INTERFACE_GRADIENT_DATA },
        { "interOverset", ONEFLOW::INTERFACE_OVERSET_DATA  }
    };

    std::string rootString =
        Prj::GetSystemFileName(
            basicString + "/alloc/" );

    OStream & logger = OStream::Instance();

    for ( const InterfaceFileSpec & spec : interfaceFileSpecs )
    {
        logger.ClearAll();
        logger << rootString << spec.name << ".txt";

        BoolIO boolIO;

        boolIO.ReadFile(
            logger.str() );

        AddInterfaceFieldNames(
            solverType,
            spec.fieldType,
            boolIO.GetNameValuePair().nameList );
    }
}

void FieldAlloc::AllocateGlobalField(
    int solverType,
    const std::string & basicString )
{
    FieldFactory::AddFieldManager( solverType );

    const FieldFileSpec fieldFileSpecs[] =
    {
        { "unsteady", FieldLocation::Inner,    true  },
        { "inner",    FieldLocation::Inner,    false },
        { "face",     FieldLocation::Face,     false },
        { "bc",       FieldLocation::Boundary, false }
    };

    std::string rootString =
        Prj::GetSystemFileName(
            basicString + "/alloc/" );

    OStream & logger = OStream::Instance();

    for ( const FieldFileSpec & spec : fieldFileSpecs )
    {
        logger.ClearAll();
        logger << rootString << spec.name << ".txt";

        ReadSuperPara readSuperPara( solverType );

        readSuperPara.Register(
            logger.str(),
            spec.location,
            spec.isUnsteady );
    }

    FieldAlloc::AllocateAllKindsOfInterfaceField(
        solverType );
}

void FieldAlloc::AllocateAllKindsOfInterfaceField( int solverType )
{
    FieldManager * fieldManager = FieldFactory::GetFieldManager( solverType );
    fieldManager->AllocateInnerAndBcField();
    FieldAlloc::AllocateInterfaceField( &fieldManager->iFieldProperty );
    FieldAlloc::AllocateOversetInterfaceField( &fieldManager->iFieldProperty );
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

void BoolIO::Add( const std::string & name, bool value )
{
    boolNameList.push_back( name );
    boolValueList.push_back( value );
}

const NameValuePair & BoolIO::GetNameValuePair() const
{
    return nameValuePair;
}

bool BoolIO::GetBoolValue(
    const std::string & varName ) const
{
    for ( int i = 0; i < boolNameList.size(); ++ i )
    {
        if ( varName == boolNameList[ i ] )
        {
            return boolValueList[ i ];
        }
    }

    Fatal( "Unknown boolean variable: " + varName );

    return false;
}

void BoolIO::ReadBool( TextFileParser & textFileParser )
{
    std::string varName =
        textFileParser.ReadNextWord();

    textFileParser.ReadNextWord();

    std::string var1 =
        textFileParser.ReadNextWord();

    std::string opName =
        textFileParser.ReadNextWord();

    std::string var2 =
        textFileParser.ReadNextWord();

    bool boolValue =
        CompareVar(
            var1,
            opName,
            var2 );

    this->Add( varName, boolValue );
}


void BoolIO::ReadSuperBool( TextFileParser & textFileParser )
{
    std::string varName =
        textFileParser.ReadNextWord();

    textFileParser.ReadNextWord();

    std::string var1 =
        textFileParser.ReadNextWord();

    std::string opName =
        textFileParser.ReadNextWord();

    std::string var2 =
        textFileParser.ReadNextWord();

    bool varValue1 =
        this->GetBoolValue( var1 );

    bool varValue2 =
        this->GetBoolValue( var2 );

    bool boolValue =
        CalcBoolLogic(
            varValue1,
            opName,
            varValue2 );

    this->Add( varName, boolValue );
}

void BoolIO::ReadName(
    TextFileParser & textFileParser )
{
    std::string varName =
        textFileParser.ReadNextWord();

    nameValuePair.nameList.push_back( varName );
}

void BoolIO::ReadNameValue(
    TextFileParser & textFileParser )
{
    std::string varName =
        textFileParser.ReadNextWord();

    nameValuePair.nameList.push_back( varName );

    Real varValue =
        textFileParser.ReadNextDigit< Real >();

    nameValuePair.valueList.push_back( varValue );
}

void BoolIO::ReadFile(
    const std::string & fileName )
{
    ReadBoolFile(
        *this,
        fileName,
        &BoolIO::ReadName );
}

void BoolIO::ReadValueFile(
    const std::string & fileName )
{
    ReadBoolFile(
        *this,
        fileName,
        &BoolIO::ReadNameValue );
}

ParaNameDim * ParaNameDimData::GetParaNameDim( FieldCategory category )
{
    switch ( category )
    {
    case FieldCategory::Common:
        return &comPara;

    case FieldCategory::Structured:
        return &strPara;

    case FieldCategory::Unstructured:
        return &unsPara;
    }

    return nullptr;
}

const ParaNameDim *
ParaNameDimData::GetParaNameDim( FieldCategory category ) const
{
    switch ( category )
    {
    case FieldCategory::Common:
        return &comPara;

    case FieldCategory::Structured:
        return &strPara;

    case FieldCategory::Unstructured:
        return &unsPara;
    }

    return nullptr;
}

ReadSuperPara::ReadSuperPara( int solverType )
    : solverType( solverType )
{
}

void ReadSuperPara::AddUnsteadyInnerFieldProperty()
{
    this->AddFieldProperties( FieldLocation::Inner );
    FieldManager * fieldManager = FieldFactory::GetFieldManager( this->solverType );

    UsdPara * usdPara = fieldManager->usdPara.get();

    const ParaNameDim * comPara =
        this->paraNameDimData.GetParaNameDim(
            FieldCategory::Common );

    int nEqu = comPara->nEquList[ 0 ];

    usdPara->Init(
        comPara->nameList,
        nEqu );
}

void ReadSuperPara::AddFieldProperties(
    FieldLocation location )
{
    FieldManager * fieldManager =
        FieldFactory::GetFieldManager( this->solverType );

    AddBasicFieldProperty(
        fieldManager,
        this->paraNameDimData.GetParaNameDim(
            FieldCategory::Unstructured ),
        location,
        FieldCategory::Unstructured );

    AddBasicFieldProperty(
        fieldManager,
        this->paraNameDimData.GetParaNameDim(
            FieldCategory::Structured ),
        location,
        FieldCategory::Structured );

    AddBasicFieldProperty(
        fieldManager,
        this->paraNameDimData.GetParaNameDim(
            FieldCategory::Common ),
        location,
        FieldCategory::Common );
}

void ReadSuperPara::Register(
    const std::string & fileName,
    FieldLocation location,
    bool isUnsteady )
{
    ReadFieldDefinitions(
        fileName,
        this->paraNameDimData );

    if ( isUnsteady )
    {
        this->AddUnsteadyInnerFieldProperty();
    }
    else
    {
        this->AddFieldProperties( location );
    }
}


EndNameSpace

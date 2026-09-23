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
#include "FieldBase.h"
#include "UsdPara.h"
#include "SolverInfo.h"
#include "SolverDef.h"
#include "TextFileParser.h"
#include "OStream.h"
#include "DataBase.h"
#include "RegisterUtils.h"
#include "Zone.h"
#include "Grid.h"
#include "UnsGrid.h"
#include "GridState.h"
#include "InterFace.h"

BeginNameSpace( ONEFLOW )

namespace
{
    struct FieldDefinition
    {
        std::string name;
        int nEqu;
    };

    int ResolveIntegerValue(
        const std::string & valueToken )
    {
        if ( Word::IsDigit( valueToken ) )
        {
            return StringToDigit< int >( valueToken );
        }

        return GetDataValue< int >( valueToken );
    }


     UsdFieldNames BuildUsdFieldNames(
        const ParaNameDim & paraNameDim )
    {
        UsdFieldNames fieldNames;

        fieldNames.flow.push_back(
            paraNameDim.GetName( 0 ) );

        fieldNames.flow.push_back(
            paraNameDim.GetName( 1 ) );

        fieldNames.flow.push_back(
            paraNameDim.GetName( 2 ) );

        fieldNames.residual.push_back(
            paraNameDim.GetName( 3 ) );

        fieldNames.residual.push_back(
            paraNameDim.GetName( 4 ) );

        fieldNames.residual.push_back(
            paraNameDim.GetName( 5 ) );

        fieldNames.dq.push_back(
            paraNameDim.GetName( 6 ) );

        return fieldNames;
    }

    struct FieldFileSpec
    {
        const char * name;
        FieldLocation location;
        bool initializeUsdPara;
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

        if ( typeName == "uns" )
        {
            return FieldCategory::Unstructured;
        }

        Fatal(
            "Unknown field category: " + typeName );

        return FieldCategory::Unstructured;
    }

    void ReadFieldDefinition(
        TextFileParser & textFileParser,
        ParaNameDimData & paraNameDimData )
    {
        FieldDefinition definition;

        definition.name =
            textFileParser.ReadNextWord();

        std::string equationCountToken =
            textFileParser.ReadNextWord();

        std::string categoryToken =
            textFileParser.ReadNextWord();

        definition.nEqu =
            ResolveIntegerValue( equationCountToken );

        FieldCategory category =
            ParseFieldCategory( categoryToken );

        ParaNameDim * paraNameDim =
            paraNameDimData.GetParaNameDim( category );

        paraNameDim->Add(
            definition.name,
            definition.nEqu );
    }

    void AddBasicFieldProperty(
        FieldManager * fieldManager,
        const ParaNameDim * paraNameDim,
        FieldLocation location,
        FieldCategory category )
    {
        int fieldCount = paraNameDim->Size();

        for ( int fieldIndex = 0; fieldIndex < fieldCount; ++ fieldIndex )
        {
            const std::string & fieldName =
                paraNameDim->GetName( fieldIndex );

            int nEqu =
                paraNameDim->GetNEqu( fieldIndex );

            fieldManager->AddField(
                fieldName,
                nEqu,
                category,
                location );
        }
    }

    void SetFieldValues(
        FieldManager * fieldManager,
        const NameValuePair & valuePair )
    {
        const int fieldCount =
            valuePair.Size();

        for ( int fieldIndex = 0;
            fieldIndex < fieldCount;
            ++ fieldIndex )
        {
            fieldManager->SetField(
                valuePair.GetName( fieldIndex ),
                valuePair.GetValue( fieldIndex ) );
        }
    }

    void AddInterfaceFieldNames(
        int solverType,
        int fieldType,
        const FieldNameList & fieldNameList )
    {
        VarNameSolver * varNameSolver =
            VarNameFactory::GetVarNameSolver(
                solverType,
                fieldType );

        int fieldCount =
            fieldNameList.Size();

        for ( int fieldIndex = 0;
            fieldIndex < fieldCount;
            ++ fieldIndex )
        {
            const std::string & fieldName =
                fieldNameList.GetName( fieldIndex );

            varNameSolver->AddFieldName( fieldName );
        }
    }

    bool CalcBoolLogic(
        bool leftValue,
        const std::string & operatorName,
        bool rightValue )
    {
        if ( operatorName == "&&" )
        {
            return leftValue && rightValue;
        }

        if ( operatorName == "||" )
        {
            return leftValue || rightValue;
        }

        Fatal(
            "Unknown boolean operator: "
            + operatorName );

        return false;
    }
    
    bool CompareValues(
        const std::string & leftToken,
        const std::string & operatorName,
        const std::string & rightToken )
    {
        int leftValue =
            ResolveIntegerValue( leftToken );

        int rightValue =
            ResolveIntegerValue( rightToken );

        if ( operatorName == ">" )
        {
            return leftValue > rightValue;
        }

        if ( operatorName == ">=" )
        {
            return leftValue >= rightValue;
        }

        if ( operatorName == "==" )
        {
            return leftValue == rightValue;
        }

        if ( operatorName == "<" )
        {
            return leftValue < rightValue;
        }

        if ( operatorName == "<=" )
        {
            return leftValue <= rightValue;
        }

        if ( operatorName == "!=" )
        {
            return leftValue != rightValue;
        }

        Fatal(
            "Unknown comparison operator: "
            + operatorName );

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

void FieldNameList::Add(
    const std::string & name )
{
    nameList.push_back( name );
}

int FieldNameList::Size() const
{
    return nameList.size();
}

const std::string & FieldNameList::GetName(
    int index ) const
{
    return nameList[ index ];
}

void NameValuePair::Add(
    const std::string & name,
    Real value )
{
    nameList.push_back( name );
    valueList.push_back( value );
}

int NameValuePair::Size() const
{
    return nameList.size();
}

const std::string & NameValuePair::GetName(
    int index ) const
{
    return nameList[ index ];
}

Real NameValuePair::GetValue(
    int index ) const
{
    return valueList[ index ];
}

void FieldAlloc::AllocateAllFields(
    int solverType,
    const std::string & basicString )
{
    FieldAlloc::RegisterInterfaceVar(
        solverType,
        basicString );

    FieldFactory::AddFieldManager(
        solverType );

    FieldManager * fieldManager =
        FieldFactory::GetFieldManager(
            solverType );

    if ( ! fieldManager->HasFieldDefinitions() )
    {
        FieldAlloc::RegisterFieldDefinitions(
            fieldManager,
            basicString );

        fieldManager->MarkFieldDefinitionsReady();
    }

    FieldAlloc::ValidateInterfaceVar(
        solverType,
        fieldManager );

    FieldAlloc::AllocateRuntimeFields(
        fieldManager );

    FieldAlloc::InitField(
        fieldManager,
        basicString );
}

void FieldAlloc::InitField(
    FieldManager * fieldManager,
    const std::string & basicString )
{
    std::string fileName = Prj::GetSystemFileName( basicString + "/alloc/init.txt" );
    BoolIO boolIO;
    boolIO.ReadValueFile( fileName );

    SetFieldValues(
        fieldManager,
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
            boolIO.GetFieldNameList() );
    }
}

void FieldAlloc::ValidateInterfaceVar(
    int solverType,
    FieldManager * fieldManager )
{
    const FieldProperty::Data & interfaceData =
        fieldManager->GetInterfaceFieldProperty().GetData();

    const int interfaceTypes[] =
    {
        ONEFLOW::INTERFACE_DATA,
        ONEFLOW::INTERFACE_DQ_DATA,
        ONEFLOW::INTERFACE_GRADIENT_DATA,
        ONEFLOW::INTERFACE_OVERSET_DATA
    };

    for ( int iType = 0; iType < 4; ++ iType )
    {
        VarNameSolver * varNameSolver =
            VarNameFactory::FindVarNameSolver(
                solverType,
                interfaceTypes[ iType ] );

        if ( varNameSolver == nullptr )
        {
            continue;
        }

        for ( int iField = 0;
            iField < varNameSolver->data.size();
            ++ iField )
        {
            const std::string & fieldName =
                varNameSolver->data[ iField ];

            if ( interfaceData.find( fieldName ) ==
                interfaceData.end() )
            {
                Fatal(
                    "Interface field is not allocated: "
                    + fieldName );
            }
        }
    }
}

void FieldAlloc::RegisterFieldDefinitions(
    FieldManager * fieldManager,
    const std::string & basicString )
{
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

        ReadSuperPara readSuperPara(
            fieldManager );

        readSuperPara.Register(
            logger.str(),
            spec.location,
            spec.initializeUsdPara );
    }

}

void FieldAlloc::AllocateRuntimeFields(
    FieldManager * fieldManager )
{
    FieldAlloc::AllocateGridFields(
        fieldManager );

    FieldAlloc::AllocateInterfaceField(
        &fieldManager->GetInterfaceFieldProperty() );

    FieldAlloc::AllocateOversetInterfaceField(
        &fieldManager->GetInterfaceFieldProperty() );
}

void FieldAlloc::AllocateGridFields(
    FieldManager * fieldManager )
{
    Grid * gridIn = Zone::GetGrid();

    if ( ONEFLOW::IsUnsGrid( gridIn->type ) )
    {
        UnsGrid * grid =
            ONEFLOW::UnsGridCast( gridIn );

        FieldAlloc::AllocateGridFields(
            grid,
            &fieldManager->GetFieldPropertyData(
                FieldCategory::Common ) );

        FieldAlloc::AllocateGridFields(
            grid,
            &fieldManager->GetFieldPropertyData(
                FieldCategory::Unstructured ) );
    }
}

void FieldAlloc::AllocateGridFields(
    UnsGrid * grid,
    FieldPropertyData * fieldPropertyData )
{
    FieldAlloc::AllocateInnerField(
        grid,
        fieldPropertyData );

    FieldAlloc::AllocateFaceField(
        grid,
        fieldPropertyData );

    FieldAlloc::AllocateBcField(
        grid,
        fieldPropertyData );
}

void FieldAlloc::AllocateInnerField(
    UnsGrid * grid,
    FieldPropertyData * fieldPropertyData )
{
    int nTCell = grid->nCells + grid->nBFaces;

    const FieldProperty::Data & data =
        fieldPropertyData->GetFieldProperty(
            FieldLocation::Inner ).GetData();

    for ( FieldProperty::Data::const_iterator iter = data.begin();
        iter != data.end();
        ++ iter )
    {
        int nTEqu = iter->second;

        ONEFLOW::CreateMRField(
            grid,
            nTEqu,
            nTCell,
            iter->first );

        MRField * field =
            ONEFLOW::GetFieldPointer< MRField >(
                grid,
                iter->first );

        ONEFLOW::ZeroField(
            field,
            nTEqu,
            nTCell );
    }
}

void FieldAlloc::AllocateFaceField(
    UnsGrid * grid,
    FieldPropertyData * fieldPropertyData )
{
    int nFaces = grid->nFaces;

    const FieldProperty::Data & data =
        fieldPropertyData->GetFieldProperty(
            FieldLocation::Face ).GetData();

    for ( FieldProperty::Data::const_iterator iter =
        data.begin();
        iter != data.end();
        ++ iter )
    {
        int nTEqu = iter->second;

        ONEFLOW::CreateMRField(
            grid,
            nTEqu,
            nFaces,
            iter->first );

        MRField * field =
            ONEFLOW::GetFieldPointer< MRField >(
                grid,
                iter->first );

        ONEFLOW::ZeroField(
            field,
            nTEqu,
            nFaces );
    }
}

void FieldAlloc::AllocateBcField(
    UnsGrid * grid,
    FieldPropertyData * fieldPropertyData )
{
    int nBFaces = grid->nBFaces;

    const FieldProperty::Data & data =
        fieldPropertyData->GetFieldProperty(
            FieldLocation::Boundary ).GetData();

    for ( FieldProperty::Data::const_iterator iter =
        data.begin();
        iter != data.end();
        ++ iter )
    {
        int nTEqu = iter->second;

        ONEFLOW::CreateMRField(
            grid,
            nTEqu,
            nBFaces,
            iter->first );

        MRField * field =
            ONEFLOW::GetFieldPointer< MRField >(
                grid,
                iter->first );

        ONEFLOW::ZeroField(
            field,
            nTEqu,
            nBFaces );
    }
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
        CompareValues(
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

    fieldNameList.Add( varName );
}

const FieldNameList & BoolIO::GetFieldNameList() const
{
    return fieldNameList;
}

const NameValuePair & BoolIO::GetNameValuePair() const
{
    return nameValuePair;
}

void BoolIO::ReadNameValue(
    TextFileParser & textFileParser )
{
    std::string varName =
        textFileParser.ReadNextWord();


    Real varValue =
        textFileParser.ReadNextDigit< Real >();

    nameValuePair.Add(
        varName,
        varValue );
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

void ParaNameDim::Add(
    const std::string & name,
    int nEqu )
{
    FieldEntry entry;

    entry.name = name;
    entry.nEqu = nEqu;

    this->fields.push_back( entry );
}

int ParaNameDim::Size() const
{
    return this->fields.size();
}

const std::string & ParaNameDim::GetName(
    int index ) const
{
    return this->fields[ index ].name;
}

int ParaNameDim::GetNEqu(
    int index ) const
{
    return this->fields[ index ].nEqu;
}

ParaNameDim * ParaNameDimData::GetParaNameDim(
    FieldCategory category )
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

    Fatal( "Unknown field category." );

    return nullptr;
}

const ParaNameDim *
ParaNameDimData::GetParaNameDim(
    FieldCategory category ) const
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

    Fatal( "Unknown field category." );

    return nullptr;
}

ReadSuperPara::ReadSuperPara(
    FieldManager * fieldManager )
    : fieldManager( fieldManager )
{
}

void ReadSuperPara::AddUnsteadyInnerFieldProperty()
{
    this->AddFieldProperties( FieldLocation::Inner );

    UsdPara * usdPara =
        &this->fieldManager->GetUsdPara();

    const ParaNameDim * comPara =
        this->paraNameDimData.GetParaNameDim(
            FieldCategory::Common );

    int nEqu =
        GetDataValue< int >( "nEqu" );

    UsdFieldNames fieldNames =
        BuildUsdFieldNames( *comPara );

    usdPara->Init(
        fieldNames.flow,
        fieldNames.residual,
        fieldNames.dq,
        nEqu );
}

void ReadSuperPara::AddFieldProperties(
    FieldLocation location )
{
    AddBasicFieldProperty(
        this->fieldManager,
        this->paraNameDimData.GetParaNameDim(
            FieldCategory::Unstructured ),
        location,
        FieldCategory::Unstructured );

    AddBasicFieldProperty(
        this->fieldManager,
        this->paraNameDimData.GetParaNameDim(
            FieldCategory::Structured ),
        location,
        FieldCategory::Structured );

    AddBasicFieldProperty(
        this->fieldManager,
        this->paraNameDimData.GetParaNameDim(
            FieldCategory::Common ),
        location,
        FieldCategory::Common );
}

void ReadSuperPara::Register(
    const std::string & fileName,
    FieldLocation location,
    bool initializeUsdPara )
{
    ReadFieldDefinitions(
        fileName,
        this->paraNameDimData );

    if ( initializeUsdPara )
    {
        this->AddUnsteadyInnerFieldProperty();
    }
    else
    {
        this->AddFieldProperties( location );
    }
}


EndNameSpace

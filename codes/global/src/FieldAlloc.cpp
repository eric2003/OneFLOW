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
    int ResolveIntegerValue(
        const std::string & valueToken )
    {
        if ( Word::IsDigit( valueToken ) )
        {
            return StringToDigit< int >( valueToken );
        }

        return GetDataValue< int >( valueToken );
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


    class FieldNameList
    {
    private:
        StringField nameList;

    public:
        void Add(
            const std::string & name );

        int Size() const;

        const std::string & GetName(
            int index ) const;
    };

    class NameValuePair
    {
    private:
        StringField nameList;
        RealField valueList;

    public:
        void Add(
            const std::string & name,
            Real value );

        int Size() const;

        const std::string & GetName(
            int index ) const;

        Real GetValue(
            int index ) const;
    };

    class BoolIO
    {
    private:
        StringField boolNameList;
        BoolField boolValueList;

        FieldNameList fieldNameList;
        NameValuePair nameValuePair;

    public:
        const FieldNameList & GetFieldNameList() const;

        const NameValuePair & GetNameValuePair() const;

        bool GetBoolValue(
            const std::string & varName ) const;

        void Add(
            const std::string & name,
            bool value );

        void ReadBool(
            TextFileParser & textFileParser );

        void ReadSuperBool(
            TextFileParser & textFileParser );

        void ReadName(
            TextFileParser & textFileParser );

        void ReadNameValue(
            TextFileParser & textFileParser );

        void ReadFile(
            const std::string & fileName );

        void ReadValueFile(
            const std::string & fileName );
    };

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

    struct FieldDefinition
    {
        std::string name;
        int nEqu;
        FieldApplicability applicability;
    };




    enum class FieldFileType
    {
        Standard,
        Unsteady
    };

    struct FieldFileSpec
    {
        const char * name;
        FieldLocation location;
        FieldFileType type;
    };

    struct InterfaceFileSpec
    {
        const char * name;
        int fieldType;
    };

    FieldApplicability ParseFieldApplicability(
        const std::string & typeName )
    {
        if ( typeName == "all" )
        {
            return FieldApplicability::All;
        }

        if ( typeName == "str" )
        {
            return FieldApplicability::Structured;
        }

        if ( typeName == "uns" )
        {
            return FieldApplicability::Unstructured;
        }

        Fatal(
            "Unknown field applicability: " + typeName );

        return FieldApplicability::Unstructured;
    }

    FieldDefinition ReadFieldDefinition(
        TextFileParser & textFileParser )
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

        definition.applicability =
            ParseFieldApplicability( categoryToken );

        return definition;
    }

    void AddFieldDefinition(
        FieldManager * fieldManager,
        const FieldDefinition & definition,
        FieldLocation location )
    {
        fieldManager->AddField(
            definition.name,
            definition.nEqu,
            definition.applicability,
            location );
    }

    void AddUsdFieldName(
        UsdFieldNames & fieldNames,
        const std::string & fieldName,
        const std::string & role )
    {
        if ( role == "flow" )
        {
            fieldNames.flow.push_back( fieldName );
            return;
        }

        if ( role == "residual" )
        {
            fieldNames.residual.push_back( fieldName );
            return;
        }

        Fatal(
            "Unknown unsteady field role: " + role );
    }

    void ReadFieldDefinitions(
        TextFileParser & textFileParser,
        FieldManager * fieldManager,
        FieldLocation location,
        UsdFieldNames * fieldNames )
    {
        while ( ! textFileParser.ReachTheEndOfFile() )
        {
            bool flag =
                textFileParser.ReadNextNonEmptyLine();

            if ( ! flag )
            {
                break;
            }

            std::string keyWord =
                textFileParser.ReadNextWord();

            if ( keyWord != "true" )
            {
                continue;
            }

            FieldDefinition definition =
                ReadFieldDefinition(
                    textFileParser );

            AddFieldDefinition(
                fieldManager,
                definition,
                location );

            if ( fieldNames == nullptr )
            {
                continue;
            }

            if ( textFileParser.NextWordIsEmpty() )
            {
                continue;
            }

            std::string role =
                textFileParser.ReadNextWord();

            AddUsdFieldName(
                *fieldNames,
                definition.name,
                role );
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
        FieldManager * fieldManager,
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

            varNameSolver->AddFieldName(
                fieldName );

            fieldManager->AddInterfaceField(
                fieldName );
        }
    }



    void RegisterFieldFile(
        FieldManager * fieldManager,
        const std::string & fileName,
        FieldLocation location,
        FieldFileType type )
    {
        TextFileParser textFileParser;

        // \t is the tab key
        std::string separator = " \r\n\t#$,;\"()";

        textFileParser.OpenFile(
            fileName,
            std::ios_base::in );

        textFileParser.SetDefaultSeparator(
            separator );

        if ( type == FieldFileType::Unsteady )
        {
            UsdFieldNames fieldNames;

            ReadFieldDefinitions(
                textFileParser,
                fieldManager,
                location,
                &fieldNames );

            UsdPara * usdPara =
                &fieldManager->GetUsdPara();

            int nEqu =
                GetDataValue< int >( "nEqu" );

            usdPara->Init(
                fieldNames.flow,
                fieldNames.residual,
                nEqu );
        }
        else
        {
            ReadFieldDefinitions(
                textFileParser,
                fieldManager,
                location,
                nullptr );
        }

        textFileParser.CloseFile();
    }

    void RegisterInterfaceVar(
        int solverType,
        FieldManager * fieldManager,
        const std::string & basicString )
    {
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
                fieldManager,
                spec.fieldType,
                boolIO.GetFieldNameList() );
        }
        fieldManager->MarkInterfaceDefinitionsReady();
    }

    void ValidateInterfaceVar(
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

    void RegisterFieldDefinitions(
        FieldManager * fieldManager,
        const std::string & basicString )
    {
        const FieldFileSpec fieldFileSpecs[] =
        {
            { "unsteady", FieldLocation::Inner,    FieldFileType::Unsteady },
            { "inner",    FieldLocation::Inner,    FieldFileType::Standard },
            { "face",     FieldLocation::Face,     FieldFileType::Standard },
            { "bc",       FieldLocation::Boundary, FieldFileType::Standard }
        };

        std::string rootString =
            Prj::GetSystemFileName(
                basicString + "/alloc/" );

        OStream & logger = OStream::Instance();

        for ( const FieldFileSpec & spec : fieldFileSpecs )
        {
            logger.ClearAll();
            logger << rootString << spec.name << ".txt";

            RegisterFieldFile(
                fieldManager,
                logger.str(),
                spec.location,
                spec.type );
        }

    }

    void AllocateInnerField(
        UnsGrid * grid,
        const FieldPropertyData * fieldPropertyData )
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

    void AllocateFaceField(
        UnsGrid * grid,
        const FieldPropertyData * fieldPropertyData )
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

    void AllocateBoundaryField(
        UnsGrid * grid,
        const FieldPropertyData * fieldPropertyData )
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

    void AllocateGridFields(
        UnsGrid * grid,
        const FieldPropertyData * fieldPropertyData )
    {
        AllocateInnerField(
            grid,
            fieldPropertyData );

        AllocateFaceField(
            grid,
            fieldPropertyData );

        AllocateBoundaryField(
            grid,
            fieldPropertyData );
    }

    void AllocateGridFields(
        FieldManager * fieldManager )
    {
        Grid * gridIn = Zone::GetGrid();

        if ( ONEFLOW::IsUnsGrid( gridIn->type ) )
        {
            UnsGrid * grid =
                ONEFLOW::UnsGridCast( gridIn );

            AllocateGridFields(
                grid,
                &fieldManager->GetFieldPropertyData(
                    FieldApplicability::All ) );

            AllocateGridFields(
                grid,
                &fieldManager->GetFieldPropertyData(
                    FieldApplicability::Unstructured ) );
        }
    }

    void AllocateInterfaceField( IFieldProperty * iFieldProperty )
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

    void AllocateOversetInterfaceField( IFieldProperty * iFieldProperty )
    {
    }

    void AllocateRuntimeFields(
        FieldManager * fieldManager )
    {
        AllocateGridFields(
            fieldManager );

        AllocateInterfaceField(
            &fieldManager->GetInterfaceFieldProperty() );

        AllocateOversetInterfaceField(
            &fieldManager->GetInterfaceFieldProperty() );
    }

    void InitField(
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


}

void FieldAlloc::AllocateAllFields(
    int solverType,
    const std::string & basicString )
{
    FieldFactory::AddFieldManager(
        solverType );

    FieldManager * fieldManager =
        FieldFactory::GetFieldManager(
            solverType );

    if ( ! fieldManager->HasFieldDefinitions() )
    {
        RegisterFieldDefinitions(
            fieldManager,
            basicString );

        fieldManager->MarkFieldDefinitionsReady();
    }

    if ( ! fieldManager->HasInterfaceDefinitions() )
    {
        RegisterInterfaceVar(
            solverType,
            fieldManager,
            basicString );

        ValidateInterfaceVar(
            solverType,
            fieldManager );
    }

    AllocateRuntimeFields(
        fieldManager );

    InitField(
        fieldManager,
        basicString );
}



EndNameSpace

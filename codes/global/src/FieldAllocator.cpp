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
#include "FieldAllocator.h"
#include "FieldAllocConfig.h"
#include "DataStorage.h"
#include "Prj.h"
#include "Fatal.h"
#include "FieldManager.h"
#include "FieldBase.h"
#include "UnsteadyFieldNames.h"
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
#include "FieldWrap.h"

BeginNameSpace( ONEFLOW )

namespace
{
    // Section A (alloc text parsing): FieldAllocConfig.h / FieldAllocConfig.cpp
    //
    // =========================================================================
    // Section B: definition registration into FieldManager
    //   inner/face/bc/unsteady.txt  -> Field definitions
    //   inter*.txt                  -> Interface Storage + communication names
    //   Validate: Communication Fields are a subset of Interface Storage
    //
    //   Everything in this namespace block only deals with parsed config text
    //   and FieldManager definitions. It never touches a live Grid or
    //   DataStorage. Runtime allocation lives in the second namespace block
    //   below (Section C).
    // =========================================================================

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

    FieldSpec ReadFieldSpec(
        TextFileParser & textFileParser )
    {
        FieldSpec definition;

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

    void RegisterFieldDefinition(
        FieldManager * fieldManager,
        const FieldSpec & definition,
        FieldLocation location )
    {
        fieldManager->AddField(
            definition.name,
            definition.nEqu,
            definition.applicability,
            location );
    }

    void AddUnsteadyFieldName(
        UnsteadyFieldNames & fieldNames,
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

    void SetUnsteadyFieldNames(
        FieldManager * fieldManager,
        const UnsteadyFieldNames & fieldNames )
    {
        fieldManager->SetUnsteadyFieldNames(
            fieldNames );
    }

    void ReadFieldDefinitions(
        TextFileParser & textFileParser,
        FieldManager * fieldManager,
        FieldLocation location,
        UnsteadyFieldNames * fieldNames )
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

            FieldSpec definition =
                ReadFieldSpec(
                    textFileParser );

            RegisterFieldDefinition(
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

            AddUnsteadyFieldName(
                *fieldNames,
                definition.name,
                role );
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
        int solverType,
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
            UnsteadyFieldNames fieldNames;

            ReadFieldDefinitions(
                textFileParser,
                fieldManager,
                location,
                &fieldNames );

            SetUnsteadyFieldNames(
                fieldManager,
                fieldNames );
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

            FieldConfigReader configReader;

            configReader.ReadFile(
                logger.str() );

            AddInterfaceFieldNames(
                solverType,
                fieldManager,
                spec.fieldType,
                configReader.GetFieldNameList() );
        }
        fieldManager->MarkInterfaceDefinitionsReady();
    }

    void ValidateCommunicationFields(
        int solverType,
        FieldManager * fieldManager )
    {
        const InterfaceFieldProperty & interfaceFieldProperty =
            fieldManager->GetInterfaceFieldProperty();

        const int interfaceTypes[] =
        {
            ONEFLOW::INTERFACE_DATA,
            ONEFLOW::INTERFACE_DQ_DATA,
            ONEFLOW::INTERFACE_GRADIENT_DATA,
            ONEFLOW::INTERFACE_OVERSET_DATA
        };

        for ( int interfaceType : interfaceTypes )
        {
            VarNameSolver * varNameSolver =
                VarNameFactory::FindVarNameSolver(
                    solverType,
                    interfaceType );

            if ( varNameSolver == nullptr )
            {
                continue;
            }

            for ( const std::string & fieldName : varNameSolver->data )
            {
                if ( ! interfaceFieldProperty.HasField( fieldName ) )
                {
                    // Definition-time check: communication name must appear
                    // in Interface Storage field list (before runtime alloc).
                    Fatal(
                        "Communication field is not in Interface Storage "
                        "definitions: "
                        + fieldName );
                }
            }
        }
    }

    void RegisterFieldDefinitions(
        int solverType,
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
                solverType,
                fieldManager,
                logger.str(),
                spec.location,
                spec.type );
        }
    }

} // end of Section B anonymous namespace

namespace
{
    // =========================================================================
    // Section C: runtime allocation on Grid / Interface DataStorage
    //   Must run after Section B. Idempotent create + post-create null check.
    //   init.txt constants applied last via InitField (Allocate pipeline step 5).
    //
    //   Everything in this namespace block touches a live Grid / DataStorage.
    //   It only reads FieldManager definitions that Section B already
    //   registered; it never parses alloc/*.txt itself.
    // =========================================================================

    FieldApplicability GetGridApplicability(
        int gridType )
    {
        if ( ONEFLOW::IsUnsGrid( gridType ) )
        {
            return FieldApplicability::Unstructured;
        }

        if ( ONEFLOW::IsStrGrid( gridType ) )
        {
            return FieldApplicability::Structured;
        }

        Fatal( "Unsupported grid type for field allocation" );

        return FieldApplicability::All;
    }

    // Idempotent MRField create on Grid or DataStorage:
    // lookup -> create if missing -> post-create null check -> ZeroField.
    template < typename Storage >
    void AllocateOneMRField(
        Storage * storage,
        const std::string & fieldName,
        int nEqu,
        int nSize,
        const char * failMessagePrefix )
    {
        // 1) Lookup before create: skip if already present (idempotent).
        MRField * field =
            ONEFLOW::GetFieldPointer< MRField >(
                storage,
                fieldName );

        if ( field != nullptr )
        {
            return;
        }

        // 2) Create and register.
        ONEFLOW::CreateMRField(
            storage,
            nEqu,
            nSize,
            fieldName );

        // 3) Lookup AFTER create: must be registered; do not ZeroField on null.
        field =
            ONEFLOW::GetFieldPointer< MRField >(
                storage,
                fieldName );

        if ( field == nullptr )
        {
            Fatal(
                std::string( failMessagePrefix )
                + fieldName );
        }

        // 4) Safe to zero: field is non-null.
        ONEFLOW::ZeroField(
            field,
            nEqu,
            nSize );
    }

    void AllocateFieldSet(
        UnsGrid * grid,
        const FieldDefinitionTable & fieldDefinition,
        int nSize )
    {
        const auto & data = fieldDefinition.GetData();

        for ( const auto & [ fieldName, nTEqu ] : data )
        {
            AllocateOneMRField(
                grid,
                fieldName,
                nTEqu,
                nSize,
                "Failed to create grid field: " );
        }
    }

    void AllocateInnerField(
        UnsGrid * grid,
        const FieldDefinitionSet * fieldDefinitions )
    {
        int nTCell = grid->nCells + grid->nBFaces;

        AllocateFieldSet(
            grid,
            fieldDefinitions->GetFieldDefinition(
                FieldLocation::Inner ),
            nTCell );
    }

    void AllocateFaceField(
        UnsGrid * grid,
        const FieldDefinitionSet * fieldDefinitions )
    {
        int nFaces = grid->nFaces;

        AllocateFieldSet(
            grid,
            fieldDefinitions->GetFieldDefinition(
                FieldLocation::Face ),
            nFaces );
    }

    void AllocateBoundaryField(
        UnsGrid * grid,
        const FieldDefinitionSet * fieldDefinitions )
    {
        int nBFaces = grid->nBFaces;

        AllocateFieldSet(
            grid,
            fieldDefinitions->GetFieldDefinition(
                FieldLocation::Boundary ),
            nBFaces );
    }

    void AllocateGridFields(
        UnsGrid * grid,
        const FieldDefinitionSet * fieldDefinitions )
    {
        AllocateInnerField(
            grid,
            fieldDefinitions );

        AllocateFaceField(
            grid,
            fieldDefinitions );

        AllocateBoundaryField(
            grid,
            fieldDefinitions );
    }

    void AllocateGridFields(
        FieldManager * fieldManager )
    {
        Grid * gridIn = Zone::GetGrid();


        if ( ONEFLOW::IsUnsGrid( gridIn->type ) )
        {
            UnsGrid * grid =
                ONEFLOW::UnsGridCast( gridIn );

            FieldApplicability applicability =
                GetGridApplicability( gridIn->type );

            // All fields are common to every supported grid type.
            AllocateGridFields(
                grid,
                &fieldManager->GetFieldDefinitionSet(
                    FieldApplicability::All ) );

            // Grid-specific fields are allocated in addition to the common fields.
            AllocateGridFields(
                grid,
                &fieldManager->GetFieldDefinitionSet(
                    applicability ) );
        }
    }

    void AllocateInterfaceField(
        InterfaceFieldProperty * interfaceFieldProperty,
        int nIFaces,
        DataStorage * dataStorage )
    {
        if ( nIFaces <= 0 ) return;

        const auto & data = interfaceFieldProperty->GetData();

        for ( const auto & [ fieldName, nTEqu ] : data )
        {
            AllocateOneMRField(
                dataStorage,
                fieldName,
                nTEqu,
                nIFaces,
                "Failed to create interface field: " );
        }
    }

    void AllocateInterfaceField( InterfaceFieldProperty * interfaceFieldProperty )
    {
        Grid * grid = Zone::GetGrid();

        InterFace * interFace = grid->interFace.get();

        if ( ! ONEFLOW::IsValid( interFace ) ) return;

        int nIFaces = grid->interFace->nIFaces;
        for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
        {
            AllocateInterfaceField( interfaceFieldProperty, nIFaces, &interFace->GetSendStorage( ghostId ) );
            AllocateInterfaceField( interfaceFieldProperty, nIFaces, &interFace->GetRecvStorage( ghostId ) );
        }
    }

    void AllocateOversetInterfaceField( InterfaceFieldProperty * interfaceFieldProperty )
    {
        // Reserved: overset interface storage allocation is not
        // implemented yet. Keep the call site in AllocateRuntimeFields
        // so the pipeline order stays stable when overset is wired in.
        (void) interfaceFieldProperty;
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

    // Used by InitField (pipeline step 5): apply init.txt constants onto
    // fields that AllocateRuntimeFields already created.
    void SetFieldValues(
        const NameValuePair & valuePair )
    {
        const int fieldCount =
            valuePair.Size();

        for ( int fieldIndex = 0;
            fieldIndex < fieldCount;
            ++ fieldIndex )
        {
            FieldHome::SetField(
                valuePair.GetName( fieldIndex ),
                valuePair.GetValue( fieldIndex ) );
        }
    }

    void InitField(
        const std::string & basicString )
    {
        // Constant initialization only: names/values come from
        // system/<basicString>/alloc/init.txt. Fields must already
        // exist (AllocateRuntimeFields ran before this step).
        std::string fileName = Prj::GetSystemFileName( basicString + "/alloc/init.txt" );
        FieldConfigReader configReader;
        configReader.ReadValueFile( fileName );

        SetFieldValues(
            configReader.GetNameValuePair() );
    }

} // end of Section C anonymous namespace

void FieldAllocator::Allocate(
    int solverType,
    const std::string & basicString )
{
    // Pipeline (order matters):
    // 1. Ensure FieldManager exists for solverType.
    // 2. Register field definitions once (inner/bc/inter/unsteady).
    // 3. Register interface field names once (inter*) and validate
    //    that Communication Fields are a subset of Interface Storage Fields.
    // 4. Allocate runtime storage on the current grid (and interface buffers).
    // 5. Apply constant values from init.txt onto already-allocated fields.

    FieldManager * fieldManager =
        FieldManagerRegistry::AddFieldManager(
            solverType );

    if ( ! fieldManager->HasFieldDefinitionSource() )
    {
        fieldManager->SetFieldDefinitionSource(
            basicString );
    }
    else if (
        fieldManager->GetFieldDefinitionSource() != basicString )
    {
        Fatal(
            "FieldManager configuration source mismatch for solverType: "
            + std::to_string( solverType ) );
    }

    if ( ! fieldManager->HasFieldDefinitions() )
    {
        RegisterFieldDefinitions(
            solverType,
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

        ValidateCommunicationFields(
            solverType,
            fieldManager );
    }

    AllocateRuntimeFields(
        fieldManager );

    InitField(
        basicString );
}


EndNameSpace

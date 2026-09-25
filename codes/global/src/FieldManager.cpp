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
#include "FieldManager.h"
#include "FieldBase.h"
#include "Fatal.h"
#include "UsdPara.h"
#include "Grid.h"
#include "GridState.h"
#include "DataBase.h"
#include "DataStorage.h"
#include "UnsGrid.h"
#include "Zone.h"
#include "ZoneState.h"
#include "InterFace.h"
#include "FaceTopo.h"
#include "UNsCom.h"

BeginNameSpace( ONEFLOW )

namespace
{
    void DumpFieldDefinition(
        std::ostream & output,
        const char * name,
        const FieldDefinitionTable & fieldDefinitions )
    {
        output
            << "  "
            << name
            << ":\n";

        if ( fieldDefinitions.Empty() )
        {
            output << "    <empty>\n";
            return;
        }

        fieldDefinitions.Dump( output );
    }

    void ValidateCompatibleFieldDefinition(
        const FieldDefinitionTable & fieldDefinitions,
        const std::string & fieldName,
        int nEqu )
    {
        if ( ! fieldDefinitions.HasField( fieldName ) )
        {
            return;
        }

        if ( fieldDefinitions.GetNEqu( fieldName ) != nEqu )
        {
            Fatal(
                "Conflicting field definition: "
                + fieldName );
        }
    }
}

void FieldDefinitionTable::AddField(
    const std::string & fieldName,
    int nEqu )
{
    // Policy:
    // - first registration: insert
    // - same name, same nEqu: idempotent (no-op)
    // - same name, different nEqu: configuration conflict
    auto iter = this->data.find( fieldName );

    if ( iter == this->data.end() )
    {
        this->data[ fieldName ] = nEqu;
        return;
    }

    if ( iter->second != nEqu )
    {
        Fatal(
            "Conflicting field definition: "
            + fieldName
            + " existing nEqu="
            + std::to_string( iter->second )
            + " new nEqu="
            + std::to_string( nEqu ) );
    }
}

bool FieldDefinitionTable::HasField(
    const std::string & fieldName ) const
{
    return this->data.find( fieldName ) != this->data.end();
}

int FieldDefinitionTable::GetNEqu(
    const std::string & fieldName ) const
{
    auto iter = this->data.find( fieldName );

    if ( iter == this->data.end() )
    {
        Fatal(
            "Field is not defined: "
            + fieldName );
    }

    return iter->second;
}

bool FieldDefinitionTable::Empty() const
{
    return this->data.empty();
}

const FieldDefinitionTable::Data & FieldDefinitionTable::GetData() const
{
    return this->data;
}

void FieldDefinitionTable::Dump(
    std::ostream & output ) const
{
    for ( FieldDefinitionTable::Data::const_iterator iter =
        this->data.begin();
        iter != this->data.end();
        ++ iter )
    {
        output
            << "    "
            << iter->first
            << "  nEqu="
            << iter->second
            << '\n';
    }
}

void InterfaceFieldProperty::AddField(
    const std::string & fieldName,
    int nEqu )
{
    this->fieldDefinitions.AddField(
        fieldName,
        nEqu );
}

bool InterfaceFieldProperty::HasField(
    const std::string & fieldName ) const
{
    return this->fieldDefinitions.HasField(
        fieldName );
}

int InterfaceFieldProperty::GetNEqu(
    const std::string & fieldName ) const
{
    return this->fieldDefinitions.GetNEqu(
        fieldName );
}

bool InterfaceFieldProperty::Empty() const
{
    return this->fieldDefinitions.Empty();
}

const FieldDefinitionTable::Data &
InterfaceFieldProperty::GetData() const
{
    return this->fieldDefinitions.GetData();
}

void InterfaceFieldProperty::Dump(
    std::ostream & output ) const
{
    this->fieldDefinitions.Dump( output );
}

void InterfaceFieldProperty::AllocateInterfaceField( int nIFaces, DataStorage * dataStorage )
{
    if ( nIFaces <= 0 ) return;

    const auto & data = this->GetData();
    for ( const auto & [ fieldName, nTEqu ] : data )
    {
        // 1) Lookup before create: skip if already present (idempotent).
        MRField * field =
            ONEFLOW::GetFieldPointer< MRField >(
                dataStorage,
                fieldName );

        if ( field != nullptr )
        {
            continue;
        }

        // 2) Create and register into the interface DataStorage.
        ONEFLOW::CreateMRField(
            dataStorage,
            nTEqu,
            nIFaces,
            fieldName );

        // 3) Lookup AFTER create: CreateMRField must have registered
        //    the field. A null here means create failed; do not call
        //    ZeroField on a null pointer.
        field =
            ONEFLOW::GetFieldPointer< MRField >(
                dataStorage,
                fieldName );

        if ( field == nullptr )
        {
            Fatal(
                "Failed to create interface field: "
                + fieldName );
        }

        // 4) Safe to zero: field is non-null.
        ONEFLOW::ZeroField(
            field,
            nTEqu,
            nIFaces );
    }
}

void InterfaceFieldProperty::UploadInterfaceValue()
{
    Grid * gridIn = Zone::GetGrid();

    if ( ! ONEFLOW::IsUnsGrid( gridIn->type ) )
    {
        return;
    }

    UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

    const auto & data = this->GetData();
    for ( const auto & [ fieldName, nEqu ] : data )
    {
        MRField * targetField =
            ONEFLOW::GetFieldPointer< MRField >(
                grid,
                fieldName );

        if ( targetField == nullptr )
        {
            Fatal(
                "Grid field is not allocated for interface upload: "
                + fieldName );
        }

        ONEFLOW::UploadInterfaceValue(
            grid,
            targetField,
            fieldName,
            nEqu );
    }
}

void InterfaceFieldProperty::DownloadInterfaceValue()
{
    Grid * gridIn = Zone::GetGrid();

    if ( ! ONEFLOW::IsUnsGrid( gridIn->type ) )
    {
        return;
    }

    UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

    const auto & data = this->GetData();
    for ( const auto & [ fieldName, nEqu ] : data )
    {
        MRField * targetField =
            ONEFLOW::GetFieldPointer< MRField >(
                grid,
                fieldName );

        if ( targetField == nullptr )
        {
            Fatal(
                "Grid field is not allocated for interface download: "
                + fieldName );
        }

        ONEFLOW::DownloadInterfaceValue(
            grid,
            targetField,
            fieldName,
            nEqu );
    }
}

void InterfaceFieldProperty::UploadOversetInterfaceValue()
{
    Grid * gridIn = Zone::GetGrid();

    if ( ONEFLOW::IsUnsGrid( gridIn->type ) )
    {
        UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

        const FieldDefinitionTable::Data & data = this->GetData();
        for ( FieldDefinitionTable::Data::const_iterator iter = data.begin(); iter != data.end(); ++ iter )
        {
            int nEqu = iter->second;

            MRField * targetField = ONEFLOW::GetFieldPointer< MRField >( grid, iter->first );

            ONEFLOW::UploadOversetValue( grid, targetField, iter->first,  nEqu );
        }
    }
}

void InterfaceFieldProperty::DownloadOversetInterfaceValue()
{
    Grid * gridIn = Zone::GetGrid();

    if ( ONEFLOW::IsUnsGrid( gridIn->type ) )
    {
        UnsGrid * grid = ONEFLOW::UnsGridCast( gridIn );

        const FieldDefinitionTable::Data & data = this->GetData();
        for ( FieldDefinitionTable::Data::const_iterator iter = data.begin(); iter != data.end(); ++ iter )
        {
            int nEqu = iter->second;

            MRField * targetField = ONEFLOW::GetFieldPointer< MRField >( grid, iter->first );

            ONEFLOW::DownloadOversetValue( grid, targetField, iter->first, nEqu );
        }
    }
}

FieldDefinitionTable & FieldDefinitionSet::GetFieldDefinition(
    FieldLocation location )
{
    switch ( location )
    {
    case FieldLocation::Inner:
        return innerField;

    case FieldLocation::Face:
        return faceField;

    case FieldLocation::Boundary:
        return bcField;
    }

    Fatal( "Invalid field location" );
    return innerField;
}

const FieldDefinitionTable & FieldDefinitionSet::GetFieldDefinition(
    FieldLocation location ) const
{
    switch ( location )
    {
    case FieldLocation::Inner:
        return innerField;

    case FieldLocation::Face:
        return faceField;

    case FieldLocation::Boundary:
        return bcField;
    }

    Fatal( "Invalid field location" );
    return innerField;
}

FieldManager::FieldManager()
    : fieldDefinitionsReady( false )
    , interfaceDefinitionsReady( false )
{
    usdPara =
        std::make_unique< UsdPara >();
}

FieldManager::~FieldManager() = default;

bool FieldManager::HasFieldDefinitions() const
{
    return this->fieldDefinitionsReady;
}

void FieldManager::MarkFieldDefinitionsReady()
{
    this->fieldDefinitionsReady = true;
}

bool FieldManager::HasInterfaceDefinitions() const
{
    return this->interfaceDefinitionsReady;
}

void FieldManager::MarkInterfaceDefinitionsReady()
{
    this->interfaceDefinitionsReady = true;
}

InterfaceFieldProperty & FieldManager::GetInterfaceFieldProperty()
{
    return this->interfaceFieldProperty;
}

const InterfaceFieldProperty & FieldManager::GetInterfaceFieldProperty() const
{
    return this->interfaceFieldProperty;
}

FieldDefinitionSet & FieldManager::GetFieldDefinitionSet(
    FieldApplicability applicability )
{
    switch ( applicability )
    {
    case FieldApplicability::All:
        return allFields;

    case FieldApplicability::Structured:
        return structuredFields;

    case FieldApplicability::Unstructured:
        return unstructuredFields;
    }

    Fatal( "Invalid field category" );
    return allFields;
}

const FieldDefinitionSet & FieldManager::GetFieldDefinitionSet(
    FieldApplicability applicability ) const
{
    switch ( applicability )
    {
    case FieldApplicability::All:
        return allFields;

    case FieldApplicability::Structured:
        return structuredFields;

    case FieldApplicability::Unstructured:
        return unstructuredFields;
    }

    Fatal( "Invalid field category" );
    return allFields;
}

UsdPara & FieldManager::GetUsdPara()
{
    return *this->usdPara;
}

const UsdPara & FieldManager::GetUsdPara() const
{
    return *this->usdPara;
}

void FieldManager::DumpFieldEnvironment(
    std::ostream & output ) const
{
    output
        << "========== Field Environment ==========\n\n";

    output
        << "[All]\n";

    DumpFieldDefinition(
        output,
        "Inner",
        this->allFields.GetFieldDefinition(
            FieldLocation::Inner ) );

    DumpFieldDefinition(
        output,
        "Face",
        this->allFields.GetFieldDefinition(
            FieldLocation::Face ) );

    DumpFieldDefinition(
        output,
        "Boundary",
        this->allFields.GetFieldDefinition(
            FieldLocation::Boundary ) );

    output
        << "\n[Structured]\n";

    DumpFieldDefinition(
        output,
        "Inner",
        this->structuredFields.GetFieldDefinition(
            FieldLocation::Inner ) );

    DumpFieldDefinition(
        output,
        "Face",
        this->structuredFields.GetFieldDefinition(
            FieldLocation::Face ) );

    DumpFieldDefinition(
        output,
        "Boundary",
        this->structuredFields.GetFieldDefinition(
            FieldLocation::Boundary ) );

    output
        << "\n[Unstructured]\n";

    DumpFieldDefinition(
        output,
        "Inner",
        this->unstructuredFields.GetFieldDefinition(
            FieldLocation::Inner ) );

    DumpFieldDefinition(
        output,
        "Face",
        this->unstructuredFields.GetFieldDefinition(
            FieldLocation::Face ) );

    DumpFieldDefinition(
        output,
        "Boundary",
        this->unstructuredFields.GetFieldDefinition(
            FieldLocation::Boundary ) );

    output
        << "\n[Interface Storage]\n";

    if ( this->interfaceFieldProperty.Empty() )
    {
        output << "    <empty>\n";
    }
    else
    {
        this->interfaceFieldProperty.Dump( output );
    }

    output
        << "\n========================================\n";
}

void FieldManager::AddField(
    const std::string & fieldName,
    int nEqu,
    FieldApplicability applicability,
    FieldLocation location )
{
    const FieldApplicability applicabilityList[] =
    {
        FieldApplicability::All,
        FieldApplicability::Structured,
        FieldApplicability::Unstructured
    };

    for ( FieldApplicability otherApplicability : applicabilityList )
    {
        // The current applicability is checked by AddField() below.
        // Here we only validate definitions from the other applicability
        // categories to keep the same field name and location consistent.
        if ( otherApplicability == applicability )
        {
            continue;
        }

        ValidateCompatibleFieldDefinition(
            this->GetFieldDefinitionSet(
                otherApplicability ).GetFieldDefinition(
                    location ),
            fieldName,
            nEqu );
    }

    FieldDefinitionTable & fieldDefinition =
        this->GetFieldDefinitionSet(
            applicability ).GetFieldDefinition(
                location );

    fieldDefinition.AddField(
        fieldName,
        nEqu );
}

void FieldManager::AddInterfaceField(
    const std::string & fieldName )
{
    int nEqu = 0;

    if ( ! this->FindInnerFieldDefinition(
        fieldName,
        nEqu ) )
    {
        Fatal(
            "Interface field is not defined in the field definitions: "
            + fieldName );
    }

    this->interfaceFieldProperty.AddField(
        fieldName,
        nEqu );
}

bool FieldManager::FindInnerFieldDefinition(
    const std::string & fieldName,
    int & nEqu ) const
{
    bool found = false;

    const FieldDefinitionSet * dataList[] =
    {
        &this->allFields,
        &this->structuredFields,
        &this->unstructuredFields
    };

    for ( const FieldDefinitionSet * fieldDefinitions : dataList )
    {
        const FieldDefinitionTable & fieldDefinition =
            fieldDefinitions->GetFieldDefinition(
                FieldLocation::Inner );

        if ( ! fieldDefinition.HasField( fieldName ) )
        {
            continue;
        }

        int fieldNEqu =
            fieldDefinition.GetNEqu( fieldName );

        if ( ! found )
        {
            nEqu = fieldNEqu;
            found = true;
        }
        else if ( nEqu != fieldNEqu )
        {
            Fatal(
                "Conflicting field definition: "
                + fieldName );
        }
    }

    return found;
}
std::map< int, std::unique_ptr< FieldManager > > FieldManagerRegistry::data;

void FieldManagerRegistry::AddFieldManager( int solverType )
{
    auto iter = FieldManagerRegistry::data.find( solverType );

    if ( iter == FieldManagerRegistry::data.end() )
    {
        FieldManagerRegistry::data[ solverType ] =
            std::make_unique< FieldManager >();
    }
}

FieldManager * FieldManagerRegistry::GetFieldManager( int solverType )
{
    auto iter = FieldManagerRegistry::data.find( solverType );

    if ( iter == FieldManagerRegistry::data.end() )
    {
        return nullptr;
    }

    return iter->second.get();
}

void FieldManagerRegistry::FreeFieldManager()
{
    FieldManagerRegistry::data.clear();
}

void UploadInterfaceValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu )
{
    InterFace * interFace = grid->interFace;
    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    int nIFaces = interFace->nIFaces;

    if ( field2D == 0 ) return;

    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        DataStorage * dataSend = interFace->dataSend[ ghostId ];

        MRField * fieldStorage =
            ONEFLOW::GetFieldPointer< MRField >( dataSend, name );

        if ( fieldStorage == nullptr )
        {
            Fatal(
                "Interface send field is not allocated: "
                + name );
        }

        for ( int iFace = 0; iFace < nIFaces; ++ iFace )
        {
            int iCell;
            grid->faceTopo->GetSId( iFace, ghostId + 1, iCell );

            for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
            {
                ( * fieldStorage )[ iEqu ][ iFace ] =
                    ( * field2D )[ iEqu ][ iCell ];
            }
        }
    }
}

void DownloadInterfaceValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu )
{
    InterFace * interFace = grid->interFace;
    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    if ( field2D == 0 ) return;

    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        DataStorage * dataRecv = interFace->dataRecv[ ghostId ];

        MRField * fieldStorage =
            ONEFLOW::GetFieldPointer< MRField >( dataRecv, name );

        if ( fieldStorage == nullptr )
        {
            Fatal(
                "Interface recv field is not allocated: "
                + name );
        }

        int nIFaces = interFace->nIFaces;
        for ( int iFace = 0; iFace < nIFaces; ++ iFace )
        {
            int iCell;
            grid->faceTopo->GetTId( iFace, ghostId + 1, iCell );

            for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
            {
                ( * field2D )[ iEqu ][ iCell ] =
                    ( * fieldStorage )[ iEqu ][ iFace ];
            }
        }
    }
}

void UploadOversetValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu )
{
    // Reserved: overset interface transfer is not implemented yet.
    (void) grid;
    (void) field2D;
    (void) name;
    (void) nEqu;
}

void DownloadOversetValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu )
{
    // Reserved: overset interface transfer is not implemented yet.
    (void) grid;
    (void) field2D;
    (void) name;
    (void) nEqu;
}

EndNameSpace

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
#include "FieldWrap.h"
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

        if ( fieldDefinitions.GetData().empty() )
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
        const FieldDefinitionTable::Data & data =
            fieldDefinitions.GetData();

        FieldDefinitionTable::Data::const_iterator iter =
            data.find( fieldName );

        if ( iter == data.end() )
        {
            return;
        }

        if ( iter->second != nEqu )
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
    FieldDefinitionTable::Data::iterator iter =
        this->data.find( fieldName );

    if ( iter == this->data.end() )
    {
        this->data[ fieldName ] = nEqu;
        return;
    }

    if ( iter->second != nEqu )
    {
        Fatal(
            "Conflicting field definition: "
            + fieldName );
    }
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

    const FieldDefinitionTable::Data & data = this->GetData();
    for ( FieldDefinitionTable::Data::const_iterator iter = data.begin(); iter != data.end(); ++ iter )
    {
        int nTEqu = iter->second;

        ONEFLOW::CreateMRField( dataStorage, nTEqu, nIFaces, iter->first );

        MRField * field = ONEFLOW::GetFieldPointer< MRField >( dataStorage, iter->first );
        ONEFLOW::ZeroField( field, nTEqu, nIFaces );
    }
}

void InterfaceFieldProperty::UploadInterfaceValue()
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
            ONEFLOW::UploadInterfaceValue( grid, targetField, iter->first,  nEqu );
        }
    }
}

void InterfaceFieldProperty::DownloadInterfaceValue()
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

            ONEFLOW::DownloadInterfaceValue( grid, targetField, iter->first,  nEqu );
        }
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

    if ( this->interfaceFieldProperty.GetData().empty() )
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

void FieldManager::SetField( const std::string & fieldName, Real value )
{
    FieldHome::SetField( fieldName, value );
}

void FieldManager::AddField(
    const std::string & fieldName,
    int nEqu,
    FieldApplicability applicability,
    FieldLocation location )
{
    if ( applicability == FieldApplicability::All )
    {
        ValidateCompatibleFieldDefinition(
            this->GetFieldDefinitionSet(
                FieldApplicability::Structured ).GetFieldDefinition(
                    location ),
            fieldName,
            nEqu );

        ValidateCompatibleFieldDefinition(
            this->GetFieldDefinitionSet(
                FieldApplicability::Unstructured ).GetFieldDefinition(
                    location ),
            fieldName,
            nEqu );
    }
    else
    {
        ValidateCompatibleFieldDefinition(
            this->GetFieldDefinitionSet(
                FieldApplicability::All ).GetFieldDefinition(
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
        const FieldDefinitionTable::Data & data =
            fieldDefinitions->GetFieldDefinition(
                FieldLocation::Inner ).GetData();

        FieldDefinitionTable::Data::const_iterator iter =
            data.find( fieldName );

        if ( iter == data.end() )
        {
            continue;
        }

        if ( ! found )
        {
            nEqu = iter->second;
            found = true;
            continue;
        }

        if ( nEqu != iter->second )
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
    std::map< int, std::unique_ptr< FieldManager > >::iterator iter =
        FieldManagerRegistry::data.find( solverType );

    if ( iter == FieldManagerRegistry::data.end() )
    {
        FieldManagerRegistry::data[ solverType ] =
            std::make_unique< FieldManager >();
    }
}

FieldManager * FieldManagerRegistry::GetFieldManager( int solverType )
{
    std::map< int, std::unique_ptr< FieldManager > >::iterator iter =
        FieldManagerRegistry::data.find( solverType );

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

        MRField * fieldStorage = ONEFLOW::GetFieldPointer< MRField >( dataSend, name );

        for ( int iFace = 0; iFace < nIFaces; ++ iFace )
        {
            int iCell;
            grid->faceTopo->GetSId( iFace, ghostId + 1, iCell );

            for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
            {
                ( * fieldStorage )[ iEqu ][ iFace ] = ( * field2D )[ iEqu ][ iCell ];
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

        MRField * fieldStorage = ONEFLOW::GetFieldPointer< MRField >( dataRecv, name );

        int nIFaces = interFace->nIFaces;
        for ( int iFace = 0; iFace < nIFaces; ++ iFace )
        {
            int iCell;
            grid->faceTopo->GetTId( iFace, ghostId + 1, iCell );

            for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
            {
                ( * field2D )[ iEqu ][ iCell ] = ( * fieldStorage )[ iEqu ][ iFace ];
            }
        }
    }
}

void DownloadInterfaceValue_TEST( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu )
{
    InterFace * interFace = grid->interFace;
    if ( ! ONEFLOW::IsValid( interFace ) ) return;

    if ( field2D == 0 ) return;

    for ( int ghostId = MAX_GHOST_LEVELS - 1; ghostId >= 0; -- ghostId )
    {
        DataStorage * dataRecv = interFace->dataRecv[ ghostId ];

        MRField * fieldStorage = ONEFLOW::GetFieldPointer< MRField >( dataRecv, name );

        int nIFaces = interFace->nIFaces;
        for ( int iFace = 0; iFace < nIFaces; ++ iFace )
        {
            int iCell;
            grid->faceTopo->GetTId( iFace, ghostId + 1, iCell );

            int iBFace = grid->interFace->i2b[ iFace ];
            int tId = grid->faceTopo->rCells[ iBFace ];

            for ( int iEqu = 0; iEqu < nEqu; ++ iEqu )
            {
                ( * field2D )[ iEqu ][ iCell ] = ( * fieldStorage )[ iEqu ][ iFace ];
            }
        }
    }
}

void UploadOversetValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu )
{
}


void DownloadOversetValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu )
{
}

EndNameSpace

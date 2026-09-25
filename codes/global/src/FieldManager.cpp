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
#include "Fatal.h"
#include "UsdPara.h"

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

EndNameSpace
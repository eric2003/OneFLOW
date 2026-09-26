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
#pragma once
#include "NamespaceMacros.h"
#include "FieldCategory.h"
#include "FieldDefinitionTable.h"
#include "InterfaceFieldProperty.h"
#include <map>
#include <memory>
#include <ostream>
#include <string>

BeginNameSpace( ONEFLOW )

class UnsGrid;

// FieldManager owns the field *definitions* for one solverType: which
// fields exist (inner/face/boundary, all/structured/unstructured) and
// which of them participate in Interface Storage communication.
// Runtime allocation on a live Grid/DataStorage is done by
// FieldAllocator (FieldAllocator.h/.cpp), not by FieldManager itself.
class FieldManager
{
public:
    FieldManager();
    ~FieldManager();
public:
    bool HasFieldDefinitions() const;
    void MarkFieldDefinitionsReady();

    bool HasInterfaceDefinitions() const;
    void MarkInterfaceDefinitionsReady();

public:
    InterfaceFieldProperty & GetInterfaceFieldProperty();

    const InterfaceFieldProperty & GetInterfaceFieldProperty() const;

    FieldDefinitionSet & GetFieldDefinitionSet(
        FieldApplicability applicability );

    const FieldDefinitionSet & GetFieldDefinitionSet(
        FieldApplicability applicability ) const;

    void DumpFieldEnvironment(
        std::ostream & output ) const;

    void AddInterfaceField(
        const std::string & fieldName );

    void AddField(
        const std::string & fieldName,
        int nEqu,
        FieldApplicability applicability,
        FieldLocation location );

private:
    bool FindInnerFieldDefinition(
        const std::string & fieldName,
        int & nEqu ) const;

    FieldDefinitionSet allFields;
    FieldDefinitionSet structuredFields;
    FieldDefinitionSet unstructuredFields;
    InterfaceFieldProperty interfaceFieldProperty;

    bool fieldDefinitionsReady;
    bool interfaceDefinitionsReady;
};

// Per-solverType FieldManager registry.
// AddFieldManager: create empty manager if missing (idempotent).
// GetFieldManager: optional lookup (nullptr if not registered).
// Callers that require a manager after Add must null-check or Fatal.
class FieldManagerRegistry
{
public:
    static FieldManager * AddFieldManager(
        int solverType );
    static FieldManager * GetFieldManager( int solverType );
    static void FreeFieldManager();

private:
    static std::map< int, std::unique_ptr< FieldManager > > data;
};

EndNameSpace
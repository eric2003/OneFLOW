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
#include "HXArray.h"
#include "FieldCategory.h"
#include <map>
#include <memory>
#include <ostream>
#include <string>

BeginNameSpace( ONEFLOW )

class FieldDefinitionTable
{
public:
    using Data = std::map< std::string, int >;

public:
    void AddField(
        const std::string & fieldName,
        int nEqu );

    bool HasField(
        const std::string & fieldName ) const;

    int GetNEqu(
        const std::string & fieldName ) const;

    bool Empty() const;

    const Data & GetData() const;

    void Dump(
        std::ostream & output ) const;

private:
    Data data;
};

class DataStorage;

class InterfaceFieldProperty
{
public:
    void AddField(
        const std::string & fieldName,
        int nEqu );

    bool HasField(
        const std::string & fieldName ) const;

    int GetNEqu(
        const std::string & fieldName ) const;

    bool Empty() const;

    const FieldDefinitionTable::Data & GetData() const;

    void Dump(
        std::ostream & output ) const;

    void AllocateInterfaceField(
        int nIFaces,
        DataStorage * dataStorage );

    void UploadInterfaceValue();
    void DownloadInterfaceValue();
    void UploadOversetInterfaceValue();
    void DownloadOversetInterfaceValue();

private:
    FieldDefinitionTable fieldDefinitions;
};

class FieldDefinitionSet
{
public:
    FieldDefinitionTable & GetFieldDefinition(
        FieldLocation location );

    const FieldDefinitionTable & GetFieldDefinition(
        FieldLocation location ) const;

private:
    FieldDefinitionTable bcField;
    FieldDefinitionTable faceField;
    FieldDefinitionTable innerField;
};

class UsdPara;
class UnsGrid;

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

    UsdPara & GetUsdPara();
    const UsdPara & GetUsdPara() const;

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
    std::unique_ptr< UsdPara > usdPara;

    bool fieldDefinitionsReady;
    bool interfaceDefinitionsReady;
};

class FieldManagerRegistry
{
public:
    static void AddFieldManager( int solverType );
    static FieldManager * GetFieldManager( int solverType );
    static void FreeFieldManager();

private:
    static std::map< int, std::unique_ptr< FieldManager > > data;
};

void UploadInterfaceValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu );
void DownloadInterfaceValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu );
void UploadOversetValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu );
void DownloadOversetValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu );

void DownloadInterfaceValue_TEST( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu );

EndNameSpace

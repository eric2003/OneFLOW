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

BeginNameSpace( ONEFLOW )

class FieldProperty
{
public:
    using Data = std::map< std::string, int >;

public:
    void AddField( const std::string & fieldName, int nEqu );

    const Data & GetData() const;

    void Dump( std::ostream & output ) const;

private:
    Data data;
};

class DataStorage;

class IFieldProperty : public FieldProperty
{
public:
    void AllocateInterfaceField( int nIFaces, DataStorage * dataStorage );
    void UploadInterfaceValue();
    void DownloadInterfaceValue();
    void UploadOversetInterfaceValue();
    void DownloadOversetInterfaceValue();
};

class FieldPropertyData
{
public:
    FieldProperty & GetFieldProperty(
        FieldLocation location );

    const FieldProperty & GetFieldProperty(
        FieldLocation location ) const;

private:
    FieldProperty bcField;
    FieldProperty faceField;
    FieldProperty innerField;
};

class FieldManager;
class UsdPara;
class FieldPropertyData;
class UnsGrid;

class FieldManager
{
public:
    FieldManager();
    ~FieldManager();
    bool HasFieldDefinitions() const;
    void MarkFieldDefinitionsReady();

public:
    IFieldProperty & GetInterfaceFieldProperty();

    const IFieldProperty & GetInterfaceFieldProperty() const;

    FieldPropertyData & GetFieldPropertyData(
        FieldCategory category );

    const FieldPropertyData & GetFieldPropertyData(
        FieldCategory category ) const;

    UsdPara & GetUsdPara();
    const UsdPara & GetUsdPara() const;

    void DumpFieldEnvironment(
        std::ostream & output ) const;

    void AddField(
        const std::string & fieldName,
        int nEqu,
        FieldCategory category,
        FieldLocation location );

    void SetField(
        const std::string & fieldName,
        Real value );

    //void AllocateGridFields();

    //void AllocateGridFields(
    //    UnsGrid * grid,
    //    FieldPropertyData * fieldPropertyData );

    //void AllocateInnerField(
    //    UnsGrid * grid,
    //    FieldPropertyData * fieldPropertyData );

    //void AllocateFaceField(
    //    UnsGrid * grid,
    //    FieldPropertyData * fieldPropertyData );

    //void AllocateBcField(
    //    UnsGrid * grid,
    //    FieldPropertyData * fieldPropertyData );
private:
    FieldPropertyData commManager;
    FieldPropertyData strManager;
    FieldPropertyData unsManager;
    IFieldProperty iFieldProperty;
    std::unique_ptr< UsdPara > usdPara;
    bool fieldDefinitionsReady;
};

class FieldFactory
{
public:
    static void AddFieldManager( int solverType );
    static FieldManager * GetFieldManager( int solverType );
    static void FreeFieldManager();

private:
    static std::map< int, std::unique_ptr< FieldManager > > data;
};

class UnsGrid;
void UploadInterfaceValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu );
void DownloadInterfaceValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu );
void UploadOversetValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu );
void DownloadOversetValue( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu );

void DownloadInterfaceValue_TEST( UnsGrid * grid, MRField * field2D, const std::string & name, int nEqu );

EndNameSpace

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

BeginNameSpace( ONEFLOW )

class FieldProperty
{
public:
    using Data = std::map< std::string, int >;

public:
    void AddField( const std::string & fieldName, int nEqu );
    int GetNEqu( const std::string & fieldName ) const;

    const Data & GetData() const;

private:
    Data data;
};

class DataStorage;

class IFieldProperty : public FieldProperty
{
public:
    void AllocateInterfaceField( int nIFaces, DataStorage * dataStorage );
    void DeAllocateInterfaceField( DataStorage * dataStorage );
    void UploadInterfaceValue();
    void DownloadInterfaceValue();
    void UploadOversetInterfaceValue();
    void DownloadOversetInterfaceValue();
};

class GFieldProperty
{
public:
    static void AddField( const std::string & fieldName, int nEqu );
    static int GetNEqu( const std::string & fieldName );
private:
    static std::map< std::string, int > data;
};

class FieldPropertyData
{
public:
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

public:
    std::unique_ptr< IFieldProperty > iFieldProperty;
    std::unique_ptr< UsdPara > usdPara;

    FieldPropertyData commManager;
    FieldPropertyData strManager;
    FieldPropertyData unsManager;

public:
    void AddField(
        const std::string & fieldName,
        int nEqu,
        FieldCategory category,
        FieldLocation location );

    void AddFaceField( const std::string & fieldName, int nEqu );
    void AddInnerField( const std::string & fieldName, int nEqu );
    void AddBcField( const std::string & fieldName, int nEqu );

    void AddInnerField(
        const std::string & fieldName,
        int nEqu,
        FieldCategory category );

    void AddFaceField(
        const std::string & fieldName,
        int nEqu,
        FieldCategory category );

    void AddBcField(
        const std::string & fieldName,
        int nEqu,
        FieldCategory category );
public:
    void SetField( const std::string & fieldName, Real value );
    void AllocateInnerAndBcField();
    void AllocateInnerAndBcField( UnsGrid * grid, FieldPropertyData * fieldPropertyData );
    void AllocateInnerField( UnsGrid * grid, FieldPropertyData * fieldPropertyData );
    void AllocateFaceField( UnsGrid * grid, FieldPropertyData * fieldPropertyData );
    void AllocateBcField( UnsGrid * grid, FieldPropertyData * fieldPropertyData );

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

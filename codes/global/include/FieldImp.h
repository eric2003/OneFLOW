/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2020 He Xin and the OneFLOW contributors.
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
#include "HXDefine.h"
#include "HXArray.h"
#include <map>
using namespace std;

BeginNameSpace( ONEFLOW )

class FieldProperty
{
public:
    FieldProperty();
    ~FieldProperty();
public:
    std::map< string, int > data;
public:
    void AddField( const string & fieldName, int nEqu );
    int GetNEqu( const string & fileName );
};

class DataStorage;

class IFieldProperty : public FieldProperty
{
public:
    IFieldProperty();
    ~IFieldProperty();
public:
    void AllocateInterfaceField( int nIFace, DataStorage * dataStorage );
    void DeAllocateInterfaceField( DataStorage * dataStorage );
    void UploadInterfaceValue();
    void DownloadInterfaceValue();
    void UploadOversetInterfaceValue();
    void DownloadOversetInterfaceValue();
};

class GFieldProperty
{
public:
    GFieldProperty();
    ~GFieldProperty();
public:
    static std::map< string, int > data;
public:
    static void AddField( const string & fieldName, int nEqu );
    static int GetNEqu( const string & fieldName );
};

class FieldPropertyData
{
public:
    FieldPropertyData();
    ~FieldPropertyData();
public:
    FieldProperty * bcField;
    FieldProperty * faceField;
    FieldProperty * innerField;
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
    IFieldProperty * iFieldProperty;
    UsdPara  * usdPara;

    FieldPropertyData * commManager;
    FieldPropertyData * strManager;
    FieldPropertyData * unsManager;
public:
    void AddFaceField( const string & fieldName, int nEqu );
    void AddInnerField( const string & fieldName, int nEqu );
    void AddBcField( const string & fieldName, int nEqu );
    void AddInnerField( const string & fieldName, int nEqu, int type );
    void AddFaceField( const string & fieldName, int nEqu, int type );
    void AddBcField( const string & fieldName, int nEqu, int type );
public:
    void SetField( const string & fieldName, Real value );
    void AllocateInnerAndBcField();
    void AllocateInnerAndBcField( UnsGrid * grid, FieldPropertyData * fieldPropertyData );
    void AllocateInnerField( UnsGrid * grid, FieldPropertyData * fieldPropertyData );
    void AllocateFaceField( UnsGrid * grid, FieldPropertyData * fieldPropertyData );
    void AllocateBcField( UnsGrid * grid, FieldPropertyData * fieldPropertyData );

};

class FieldFactory
{
public:
    FieldFactory();
    ~FieldFactory();
public:
    static map< int, FieldManager * > * data;
public:
    static void Init();
    static void AddFieldManager( int sTid );
    static FieldManager * GetFieldManager( int sTid );
    static void FreeFieldManager();
};

class UnsGrid;
void UploadInterfaceValue( UnsGrid * grid, MRField * field2D, const string & name, int nEqu );
void DownloadInterfaceValue( UnsGrid * grid, MRField * field2D, const string & name, int nEqu );
void UploadOversetValue( UnsGrid * grid, MRField * field2D, const string & name, int nEqu );
void DownloadOversetValue( UnsGrid * grid, MRField * field2D, const string & name, int nEqu );

void DownloadInterfaceValue_TEST( UnsGrid * grid, MRField * field2D, const string & name, int nEqu );

EndNameSpace
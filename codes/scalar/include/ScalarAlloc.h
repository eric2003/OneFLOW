/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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
#include "FieldBase.h"
#include <map>


BeginNameSpace( ONEFLOW )

class ScalarFieldProperty
{
public:
    ScalarFieldProperty();
    ~ScalarFieldProperty();
public:
    std::map< std::string, int > data;
public:
    void AddField( const std::string & fieldName, int nEqu );
    int GetNEqu( const std::string & fileName );
};

class DataStorage;

class ScalarFieldAlloc : public ScalarFieldProperty
{
public:
    ScalarFieldAlloc();
    ~ScalarFieldAlloc();
public:
    template < typename TStorage >
    void AllocateField( TStorage * storage, int nElements )
    {
        if ( nElements <= 0 ) return;

        for ( std::map< std::string, int >::iterator iter = this->data.begin(); iter != this->data.end(); ++ iter )
        {
            int nTEqu = iter->second;

            ONEFLOW::CreateMRField( storage, nTEqu, nElements, iter->first );

            MRField * field = ONEFLOW::GetFieldPointer< MRField >( storage, iter->first );
            ONEFLOW::ZeroField( field, nTEqu, nElements );
        }
    }

    template < typename TStorage >
    void DeAllocateField( TStorage * storage )
    {
        for ( std::map< std::string, int >::iterator iter = this->data.begin(); iter != this->data.end(); ++ iter )
        {
        }
    }
};


class ScalarFieldManager
{
public:
    ScalarFieldManager();
    ~ScalarFieldManager();
public:
    ScalarFieldAlloc * interfaceAlloc;
    ScalarFieldAlloc * inner;
    ScalarFieldAlloc * faceField;
public:
    void Init();
    void AllocateAllFields();
    void AllocateInterfaceField();
    void AllocateFaceField();
    void AllocateInnnerField();
    void UploadInterfaceField();
    void DownloadInterfaceField();
};

class ScalarGrid;

void ScalarUploadInterfaceValue( ScalarGrid * grid, const std::string & name );
void ScalarDownloadInterfaceValue( ScalarGrid * grid, const std::string & name );

EndNameSpace

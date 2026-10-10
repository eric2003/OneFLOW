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
#include "HXDefine.h"
#include "HXArray.h"
#include <map>


BeginNameSpace( ONEFLOW )

class DataStorage;

// Static name -> equation-count table (process-wide).
class GFieldDim
{
public:
    GFieldDim() = default;
    ~GFieldDim() = default;
public:
    static std::map< std::string, int > data;
public:
    static void AddField( const std::string & fieldName, int nEqu );
    static int GetNEqu( const std::string & fieldName );
};

// Non-owning view of MRField pointers plus parallel nEqu list.
// Callers (DataStorage / grid) keep field lifetime.
class ScalarFieldRecord
{
public:
    ScalarFieldRecord() = default;
    ~ScalarFieldRecord() = default;
public:
    void AddField( MRField * field, int nEqu );
    MRField * GetField( int id );
    const MRField * GetField( int id ) const;
    int GetNumberOfFields() const
    {
        return static_cast< int >( this->fields.size() );
    }
    int GetNEqu( int id ) const;
public:
    void AddFieldRecord( DataStorage * dataStorage, StringField & fieldNameList );
private:
    IntField nEquList;
    HXVector< MRField * > fields;
};


EndNameSpace

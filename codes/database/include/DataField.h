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
#include <unordered_map>
#include <string>

BeginNameSpace( ONEFLOW )

class PointerWrap;

class FieldEntry
{
public:
    FieldEntry();
    FieldEntry( const std::string & name, PointerWrap * data );
    ~FieldEntry();
public:
    std::string   name;
    PointerWrap * data;
public:
    std::string & GetName() { return name; }
    PointerWrap * GetPointerWrap() { return data; }
};

class DataField
{
public:
    // Use unordered_map for O(1) average lookup
    using DataMap = std::unordered_map<std::string, FieldEntry*>;
public:
    DataField();
    ~DataField();
protected:
    DataMap * dataMap;
public:
    void UpdateFieldEntry( FieldEntry * fieldEntry );
    FieldEntry * GetFieldEntry( const std::string & name );
    void DeleteFieldEntry( const std::string & name );

    DataMap * GetDataMap() { return dataMap; }
};

EndNameSpace

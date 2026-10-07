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
#include <memory>
#include <unordered_map>
#include <string>

BeginNameSpace( ONEFLOW )

class PointerWrap;

class FieldEntry
{
public:
    FieldEntry();
    FieldEntry( const std::string & name, std::unique_ptr<PointerWrap> data );
    ~FieldEntry();
public:
    std::string name;
    std::unique_ptr<PointerWrap> data;
public:
    std::string & GetName() { return name; }
    const std::string & GetName() const { return name; }
    PointerWrap * GetPointerWrap() { return data.get(); }
    const PointerWrap * GetPointerWrap() const { return data.get(); }
};

class DataField
{
public:
    using DataMap = std::unordered_map<std::string, std::unique_ptr<FieldEntry>>;
public:
    DataField();
    ~DataField();
protected:
    DataMap dataMap;
public:
    // Takes ownership of fieldEntry.
    void UpdateFieldEntry( std::unique_ptr<FieldEntry> fieldEntry );
    FieldEntry * GetFieldEntry( const std::string & name );
    const FieldEntry * GetFieldEntry( const std::string & name ) const;
    void DeleteFieldEntry( const std::string & name );
    void Clear();

    DataMap * GetDataMap() { return &dataMap; }
    const DataMap * GetDataMap() const { return &dataMap; }
};

EndNameSpace

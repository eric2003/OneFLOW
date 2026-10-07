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
#include <fstream>

BeginNameSpace( ONEFLOW )

class DataObject;

class DataEntry
{
public:
    DataEntry( const std::string & name, int type, int size, std::unique_ptr<DataObject> data );
    ~DataEntry();
private:
    const std::string name;
    const int type;
    const int size;
    std::unique_ptr<DataObject> data;
public:
    const std::string & GetName() const { return name; }
    int GetType() const { return type; }
    int GetSize() const { return size; }
    DataObject * GetDataObject() { return data.get(); }
    const DataObject * GetDataObject() const { return data.get(); }

    void Copy( const DataEntry & inputData );
    void Dump( std::fstream & file ) const;
};

class DataPara
{
public:
    // Use unordered_map for O(1) average lookup
    using DataMap = std::unordered_map< std::string, std::unique_ptr<DataEntry> >;
public:
    DataPara();
    ~DataPara();
protected:
    DataMap dataMap;
public:
    // Takes ownership of data.
    void SetDataEntry( std::unique_ptr<DataEntry> data );
    DataEntry * GetDataPointer( const std::string & name );
    const DataEntry * GetDataPointer( const std::string & name ) const;
    void RemoveDataEntry( const std::string & name );

    // Release all case-local parameter entries while keeping the database alive.
    void Clear();

    const DataMap & GetDataMap() const { return dataMap; }

    void DumpData( std::fstream & file ) const;
};

EndNameSpace

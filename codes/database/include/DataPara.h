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
#include <fstream>

BeginNameSpace( ONEFLOW )

class DataObject;

class DataV
{
public:
    DataV();
    DataV( const std::string & name, int type, int size, DataObject * data );
    ~DataV();
public:
    std::string  name;
    int          type;
    int          size;
    DataObject * data;
public:
    void Copy( DataV * inputData );
    void Dump( std::fstream & file );
};

class DataPara
{
public:
    // Use unordered_map for O(1) average lookup
    typedef std::unordered_map< std::string, DataV * > DataMap;
public:
    DataPara();
    ~DataPara();
protected:
    DataMap * dataMap;
public:
    void UpdateDataPointer( DataV * data );
    DataV * GetDataPointer( const std::string & name );
    void DeleteDataPointer( const std::string & name );

    // Keep old name as alias for compatibility
    DataMap * GetDataSet() { return dataMap; }
    DataMap * GetDataMap() { return dataMap; }

    void DumpData( std::fstream & file );
};

EndNameSpace

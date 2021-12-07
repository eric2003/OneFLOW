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
#include "Configure.h"
#include <set>
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
    int     type;
    int     size;
    DataObject * data;
public:
    void Copy( DataV * inputData );
    void Dump( std::fstream & file );
};

class CompareDataV
{
public:
    bool operator()( const DataV * lhs, const DataV * rhs ) const
    {
        return lhs->name < rhs->name;
    }
};

class DataPara
{
public:
    DataPara();
    ~DataPara();
public:
    typedef std::set< DataV *, CompareDataV > DataSET;
protected:
    DataSET * dataSet;
public:
    void UpdateDataPointer( DataV * data );
    DataV * GetDataPointer( const std::string & name );
    void DeleteDataPointer( const std::string & name );

    DataSET * GetDataSet() { return dataSet; }
public:
    void DumpData( std::fstream & file );
};

EndNameSpace

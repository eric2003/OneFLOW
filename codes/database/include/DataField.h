/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2023 He Xin and the OneFLOW contributors.
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

BeginNameSpace( ONEFLOW )

class PointerWrap;

class DataF
{
public:
    DataF();
    DataF( const std::string & name, PointerWrap * data );
    ~DataF();
public:
    std::string  name;
    PointerWrap * data;
public:
    std::string & GetName() { return name;  }
    PointerWrap * GetPointerWrap() { return data;  }
};

class CompareDataF
{
public:
    bool operator()( const DataF * lhs, const DataF * rhs ) const
    {
        return lhs->name < rhs->name;
    }
};

class DataField
{
public:
    typedef std::set< DataF *, CompareDataF > DataSET;
public:
    DataField();
    ~DataField();
protected:
    DataSET * dataSet;
public:
    void UpdateDataF( DataF * dataf );
    DataF * GetDataF( const std::string & name );
    void DeleteDataF( const std::string & name );

    DataSET * GetDataSet() { return dataSet; }
};

EndNameSpace

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
#include <map>
#include <string>
using namespace std;

BeginNameSpace( ONEFLOW )

class VarNameSolver
{
public:
    VarNameSolver();
    ~VarNameSolver();
public:
    StringField data;
public:
    void AddFieldName( const string & fieldName );
};

class MapIntInt;

class VarNameFactory
{
public:
    VarNameFactory();
    ~VarNameFactory();
public:
    static map< int, VarNameSolver * > * data;
    static MapIntInt * mapData;
public:
    static void Init();
    static void AddVarNameSolver( int a, int b );
    static VarNameSolver * GetVarNameSolver( int a, int b );
    static void FreeVarNameSolver();
};

class DataAB
{
public:
    DataAB(){};
    ~DataAB(){};
public:
    int a, b;
};

class CmpDataAB
{
public:
    bool operator()( const DataAB & k1, const DataAB & k2 ) const;
};

class MapIntInt
{
public:
    MapIntInt();
    ~MapIntInt();
public:
    map< DataAB, int, CmpDataAB > data;
public:
    void AddData( int a, int b );
    int  GetId( int a, int b );
};

EndNameSpace
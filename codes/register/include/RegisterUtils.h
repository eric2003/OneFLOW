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
#include <map>
#include <string>
#include <ostream>

BeginNameSpace( ONEFLOW )

// Communication field-name list for one (solverType, interface group) pair.
// Populated from system/<solver>/alloc/inter*.txt (names only; nEqu comes
// from FieldManager inner definitions). Used by interface send/recv packing.
class VarNameSolver
{
public:
    VarNameSolver();
    ~VarNameSolver();
public:
    StringField data;
public:
    void AddFieldName( const std::string & fieldName );
};

class MapIntInt;

// Registry of VarNameSolver keyed by (solverType, interface group).
// AddVarNameSolver: create empty slot (called at solver registration).
// GetVarNameSolver: required lookup (Fatal if missing).
// FindVarNameSolver: optional lookup (nullptr if missing).
class VarNameFactory
{
public:
    VarNameFactory();
    ~VarNameFactory();

public:
    static std::map< int, VarNameSolver * > * data;
    static MapIntInt * mapData;

public:
    static void Init();
    static void AddVarNameSolver( int a, int b );
    static VarNameSolver * GetVarNameSolver( int a, int b );
    static void FreeVarNameSolver();

    static VarNameSolver * FindVarNameSolver(
        int a,
        int b );

    static void Dump(
        std::ostream & output,
        int solverType );
};

// Composite key: a = solverType, b = interface group (INTERFACE_*).
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

// Maps DataAB -> dense id used as index into VarNameFactory::data.
// GetId requires a prior AddData (Fatal if key is missing).
class MapIntInt
{
public:
    MapIntInt();
    ~MapIntInt();
public:
    std::map< DataAB, int, CmpDataAB > data;
public:
    void AddData( int a, int b );
    int  GetId( int a, int b );
};

EndNameSpace

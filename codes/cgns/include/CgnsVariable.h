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
#include "HXCgns.h"
#include <string>
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

typedef double VEC_DATA;

class CgnsVector
{
public:
    CgnsVector();
    ~CgnsVector();
public:
    CGNS_ENUMT( DataType_t ) dataType;
    char name[ 33 ];
    int arrayId;
    int ndim;
    cgsize_t dims[ 12 ];
    VEC_DATA * data;
public:
    void ReadArray();
    void ReadArrayInfo( int arrayId );
    void ReadArrayContent();
    void PrintData();
public:
    void Create();
};

class CgnsZVector
{
public:
    CgnsZVector();
    ~CgnsZVector();
public:
    vector< CgnsVector * > cgnsVectorList;
    void ReadArray( int nArrays );
    void ReadArray();
};

class CgnsBase;
class CgnsUserData
{
public:
    CgnsUserData( CgnsBase * cgnsBase );
    ~CgnsUserData();
public:
    CgnsBase * cgnsBase;
    vector< CgnsZVector * > cgnsZVectorList;
public:
    void ReadUserData();
};

//class CgnsVariable
//{
//public:
//    CgnsVariable();
//    ~CgnsVariable();
//public:
//    char name[33];
//    int type;
//    int id;
//    int valid;
//    int len;
//    int datatype;
//    int dataclass;
//    int hasunits;
//    int units[5];
//    int hasconv;
//    double dataconv[2];
//    int hasexp;
//    double exponent[5];
//    CgnsVector *vd;
//};

EndNameSpace
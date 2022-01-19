/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2022 He Xin and the OneFLOW contributors.
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
#include <string>
#include <map>

BeginNameSpace( ONEFLOW )

const int HX_INT    = 1;
const int HX_FLOAT  = 2;
const int HX_DOUBLE = 3;
const int HX_REAL   = 4;
const int HX_STRING = 5;
const int HX_BOOL   = 6;

class DataBaseType
{
public:
    DataBaseType();
    ~DataBaseType();
public:
    static std::map< int, std::string > nameMap;
    static std::map< std::string, int > indexMap;
    static bool init_flag;
public:
    static void Init();
    static void AddAllItem();
    static void AddItem( const std::string &name, int index );
    static int GetIndex( const std::string & name );
    static std::string & GetName( int index );
};




EndNameSpace

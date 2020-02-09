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
#include "Configure.h"
#include <map>
#include <string>
using namespace std;

BeginNameSpace( ONEFLOW )

class ActionMapImp;

class ActionMap
{
public:
    ActionMap();
    ~ActionMap();
public:
    static ActionMapImp * imp;
public:
    static int    GetActionId( const string & name );
    static string GetActionName( int id );
    static void ReadFile( const string & fileName );
public:
    static void Init();
    static void Free();
};

class ActionMapImp
{
public:
    ActionMapImp();
    ~ActionMapImp();
public:
    map< string, int > * nameMap;
    map< int, string > * idMap;
public:
    void Register( const string & name );
    void Unregister( const string & name );
    int    GetActionId( const string & name );
    string GetActionName( int id );
    void ReadFile( const string & fileName );
};

EndNameSpace
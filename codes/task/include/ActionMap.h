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
#include <map>
#include <vector>
#include <string>



BeginNameSpace( ONEFLOW )

//class ActionMapImp;
//
//class ActionMap
//{
//public:
//    ActionMap();
//    ~ActionMap();
//public:
//    static ActionMapImp * imp;
//public:
//    static int    GetActionId( const std::string & name );
//    static std::string GetActionName( int id );
//    static void ReadFile( const std::string & fileName );
//public:
//    static void Init();
//    static void Free();
//};
//
//class ActionMapImp
//{
//public:
//    ActionMapImp();
//    ~ActionMapImp();
//public:
//    std::map< std::string, int > * nameMap;
//    std::map< int, std::string > * idMap;
//public:
//    void Register( const std::string & name );
//    void Unregister( const std::string & name );
//    int    GetActionId( const std::string & name );
//    std::string GetActionName( int id );
//    void ReadFile( const std::string & fileName );
//};


class ActionMapImp
{
public:
    // Rule of Zero: no heap-allocated members, so no need for a custom
    // constructor/destructor at all.
    ActionMapImp() = default;
    ~ActionMapImp() = default;

public:
    void Register( const std::string & name );
    void Unregister( const std::string & name );
    int GetActionId( const std::string & name ) const;
    std::string GetActionName( int id ) const;
    void ReadFile( const std::string & fileName );

    // Resets to a fresh, empty state. Used by ActionMap::Init()/Free()
    // (see ActionMap below) and directly useful for test isolation.
    void Clear();

private:
    // name -> sequential id, assigned in registration order
    std::map< std::string, int > nameToId;

    // id -> name, indexed directly. Since ids are always assigned as
    // 0, 1, 2, ... in registration order (see Register()), a vector
    // gives O(1) lookup with better cache locality than the previous
    // std::map<int, std::string>, which paid for tree lookups on keys
    // that are never actually sparse.
    std::vector< std::string > idToName;
};

class ActionMap
{
public:
    static int GetActionId( const std::string & name );
    static std::string GetActionName( int id );
    static void Register( const std::string & name );
    static void ReadFile( const std::string & fileName );

    // Kept for source compatibility with existing call sites. The
    // underlying implementation is now a self-managing function-local
    // static (see GetImp()) with no manual new/delete, so these no
    // longer allocate/deallocate anything - they just reset state to
    // empty. Safe to call any number of times, in any order.
    static void Init();
    static void Free();

private:
    static ActionMapImp & GetImp();
};

EndNameSpace

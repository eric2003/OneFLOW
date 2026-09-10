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

//
//#pragma once
//#include "HXDefine.h"
//#include <map>
//#include <string>
//
//
//BeginNameSpace( ONEFLOW )
//class HXClone;
//
//class HXRegister
//{
//public:
//    HXRegister();
//    ~HXRegister();
//public:
//    std::map< std::string, HXClone * > data;
//public:
//    void FreeAll();
//    void Register( const std::string & cmdName, const std::string & className );
//    HXClone * GetClass( const std::string & cmdName );
//};
//
//class MRegister
//{
//public:
//    MRegister();
//    ~MRegister();
//public:
//    std::vector< HXRegister * > data;
//    StringField fileNames;
//public:
//    void SetSolverFileNames( StringField & fileNames );
//public:
//    HXRegister * GetRegister( int index );
//    HXRegister * GetRegister();
//    void RegisterAll();
//private:
//    void AllocateData();
//    void Register( const std::string & fileName, HXRegister * fRegister );
//};
//
//
//class RegisterFactory
//{
//public:
//    RegisterFactory();
//    ~RegisterFactory();
//public:
//    static std::map< int, MRegister * > * data;
//public:
//    static void Init();
//    static void AddMRegister( int registerId );
//    static MRegister * GetMRegister( int registerId );
//    static void FreeMRegister();
//public:
//    static HXRegister * GetRegister( int mRegisterId, int registerId );
//};
//
//EndNameSpace

#pragma once
#include "HXDefine.h"
#include <map>
#include <memory>
#include <string>
#include <vector>

BeginNameSpace( ONEFLOW )
class HXClone;

class HXRegister
{
public:
    HXRegister() = default;
    // FIX: previously an empty destructor left every HXClone* in `data`
    // leaked, since only FreeAll() (never automatically called) freed
    // them. Switching to unique_ptr means destruction is automatic and
    // exception-safe, and this class no longer needs a hand-written
    // destructor at all (Rule of Zero).
    ~HXRegister() = default;

public:
    void Register( const std::string & cmdName, const std::string & className );
    HXClone * GetClass( const std::string & cmdName );

    // Kept for source compatibility with existing call sites that may
    // call FreeAll() explicitly; now just an alias for clearing the map,
    // which releases all owned HXClone objects via unique_ptr.
    void FreeAll();

private:
    std::map< std::string, std::unique_ptr< HXClone > > data;
};

class MRegister
{
public:
    MRegister() = default;
    // FIX: unique_ptr in `data` makes destruction automatic; no more
    // hand-written loop that only freed the HXRegister* itself while
    // silently leaking everything each HXRegister owned.
    ~MRegister() = default;

public:
    std::vector< std::unique_ptr< HXRegister > > data;
    StringField fileNames;
public:
    // Directly installs a fully-constructed HXRegister at the given
    // index, growing `data` as needed. Bypasses RegisterAll()'s
    // file-based loading entirely. Useful both for programmatically
    // wiring up fixed/built-in registrations and for tests that need
    // to inject a stub class without real config files on disk.
    void SetRegister( int index, std::unique_ptr< HXRegister > reg );

public:
    void SetSolverFileNames( StringField & fileNames );

public:
    // Returns nullptr on an out-of-range index instead of invoking
    // undefined behavior via vector::operator[]. Callers must check
    // for nullptr (this is a deliberate, minimal-footprint fix - see
    // notes below on why a broader API change is deferred).
    HXRegister * GetRegister( int index );
    HXRegister * GetRegister();
    void RegisterAll();

private:
    void AllocateData();
    void Register( const std::string & fileName, HXRegister * fRegister );
};

class RegisterFactory
{
public:
    RegisterFactory() = default;
    ~RegisterFactory() = default;

public:
    static void Init();
    static void AddMRegister( int registerId );

    // Returns nullptr if RegisterFactory::Init() was never called, or if
    // registerId was never added via AddMRegister(). Callers must check
    // for nullptr.
    static MRegister * GetMRegister( int registerId );
    static void FreeMRegister();

public:
    // Returns nullptr if either lookup fails, instead of dereferencing
    // a null MRegister* or an end() iterator.
    static HXRegister * GetRegister( int mRegisterId, int registerId );

private:
    static std::map< int, std::unique_ptr< MRegister > > & GetData();
};

EndNameSpace

/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
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


BeginNameSpace( ONEFLOW )
class HXClone;

class HXRegister
{
public:
    HXRegister();
    ~HXRegister();
public:
    std::map< std::string, HXClone * > data;
public:
    void FreeAll();
    void Register( const std::string & cmdName, const std::string & className );
    HXClone * GetClass( const std::string & cmdName );
};

class MRegister
{
public:
    MRegister();
    ~MRegister();
public:
    std::vector< HXRegister * > data;
    StringField fileNames;
public:
    void SetSolverFileNames( StringField & fileNames );
public:
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
    RegisterFactory();
    ~RegisterFactory();
public:
    static std::map< int, MRegister * > * data;
public:
    static void Init();
    static void AddMRegister( int registerId );
    static MRegister * GetMRegister( int registerId );
    static void FreeMRegister();
public:
    static HXRegister * GetRegister( int mRegisterId, int registerId );
};

EndNameSpace

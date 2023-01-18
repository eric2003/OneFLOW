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
#include "HXDefine.h"
#include <map>
#include <string>


BeginNameSpace( ONEFLOW )

class MessageMap
{
public:
    MessageMap();
    ~MessageMap();
public:
    static std::map< std::string, int > * nameMap;
    static std::map< int, std::string > * idMap;
public:
    static void Register( const std::string & msgName );
    static void Unregister( const std::string & msgName );
    static int    GetMsgId( const std::string & msgName );
    static std::string GetMsgName( int msgId );
    static void ReadFile( const std::string & fileName );
public:
    static void Init();
    static void Free();
};

EndNameSpace

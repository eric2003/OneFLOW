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


//#pragma once
//#include "NamespaceMacros.h"
//#include "HXDefine.h"
//#include <map>
//#include <string>
//
//
//BeginNameSpace( ONEFLOW )
//
//class MessageMap
//{
//public:
//    MessageMap();
//    ~MessageMap();
//public:
//    static std::map< std::string, int > * nameMap;
//    static std::map< int, std::string > * idMap;
//public:
//    static void Register( const std::string & msgName );
//    static void Unregister( const std::string & msgName );
//    static int    GetMsgId( const std::string & msgName );
//    static std::string GetMsgName( int msgId );
//    static void ReadFile( const std::string & fileName );
//public:
//    static void Init();
//    static void Free();
//};
//
//EndNameSpace


#pragma once
#include "NamespaceMacros.h"
#include "HXDefine.h"
#include <map>
#include <vector>
#include <string>

BeginNameSpace( ONEFLOW )

// Implementation class: holds the actual name<->id mapping. Instantiable
// on its own (no static state), which makes it directly unit-testable
// without going through the static MessageMap facade below.
class MessageMapImp
{
public:
    MessageMapImp() = default;
    ~MessageMapImp() = default;

public:
    void Register( const std::string & msgName );
    void Unregister( const std::string & msgName );
    int GetMsgId( const std::string & msgName ) const;
    std::string GetMsgName( int msgId ) const;
    void ReadFile( const std::string & fileName );
    void Clear();

private:
    std::map< std::string, int > nameToId;

    // ids are assigned sequentially (0, 1, 2, ...) in registration order,
    // so a vector gives O(1) lookup instead of paying for a red-black
    // tree lookup on keys that are never actually sparse.
    std::vector< std::string > idToName;
};

class MessageMap
{
public:
    static int GetMsgId( const std::string & msgName );
    static std::string GetMsgName( int msgId );
    static void Register( const std::string & msgName );
    static void ReadFile( const std::string & fileName );

    // Kept for source compatibility. No manual new/delete underneath any
    // more (see GetImp()), so these just reset state to empty and are
    // safe to call any number of times, in any order.
    static void Init();
    static void Free();

private:
    static MessageMapImp & GetImp();
};

EndNameSpace
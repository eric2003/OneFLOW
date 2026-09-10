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
#include "ActionMap.h"
#include "TextFileParser.h"

BeginNameSpace( ONEFLOW )

void ActionMapImp::Register( const std::string & actionName )
{
    // Defensive guard: an empty name should never be a valid registration.
    // This also protects against any future caller path that might pass
    // through an empty token by accident (a category of bug we already
    // found once, in ReadFile's old EOF-handling logic).
    if ( actionName.empty() )
    {
        return;
    }

    if ( this->nameToId.find( actionName ) != this->nameToId.end() )
    {
        return; // already registered; first registration wins (unchanged behavior)
    }

    int actionId = static_cast< int >( this->idToName.size() );
    this->nameToId[ actionName ] = actionId;
    this->idToName.push_back( actionName );
}

void ActionMapImp::Unregister( const std::string & actionName )
{
    // NOTE (pre-existing behavior, unchanged by this refactor): this only
    // removes the name->id mapping and does NOT remove the corresponding
    // slot from idToName (the original idMap had the same gap). A stale
    // id->name mapping therefore remains reachable via GetActionName()
    // after Unregister(). This is a known, separate issue - flagged here
    // rather than silently fixed, since no caller or test currently
    // exercises Unregister() and fixing it changes id-stability semantics
    // that deserve their own discussion and tests.
    this->nameToId.erase( actionName );
}

int ActionMapImp::GetActionId( const std::string & actionName ) const
{
    auto iter = this->nameToId.find( actionName );
    if ( iter == this->nameToId.end() )
    {
        return -1;
    }
    return iter->second;
}

std::string ActionMapImp::GetActionName( int actionId ) const
{
    if ( actionId < 0 || actionId >= static_cast< int >( this->idToName.size() ) )
    {
        return "";
    }
    return this->idToName[ actionId ];
}

void ActionMapImp::ReadFile( const std::string & fileName )
{
    std::string separator = " =\r\n\t#$,;\"";

    TextFileParser textFileParser;
    textFileParser.OpenFile( fileName, std::ios_base::in );
    textFileParser.SetDefaultSeparator( separator );

    // Driven by ReadNextMeaningfulLine()'s return value (skips blank AND
    // comment lines, and reliably signals "no more content" at EOF)
    // rather than a separate ReachTheEndOfFile() pre-check. See the fix
    // history from the previous refactor round.
    while ( textFileParser.ReadNextMeaningfulLine() )
    {
        std::string actionName = textFileParser.ReadNextWord();
        this->Register( actionName );
    }

    textFileParser.CloseFile();
}

void ActionMapImp::Clear()
{
    this->nameToId.clear();
    this->idToName.clear();
}

ActionMapImp & ActionMap::GetImp()
{
    // Meyer's singleton: constructed on first use, thread-safe since
    // C++11, destroyed automatically at program exit. This eliminates,
    // by construction, the whole class of bugs previously found in
    // Category and TaskRegister (use before Init(), dangling pointer
    // after Free(), double-free on repeated Free() calls): there is no
    // pointer and no manual new/delete here at all.
    static ActionMapImp imp;
    return imp;
}

void ActionMap::Init()
{
    // Historically allocated a fresh ActionMapImp. The implementation is
    // now a self-managing singleton, so Init() just guarantees a clean,
    // empty starting state for callers (including tests) that rely on
    // calling Init() before first use.
    ActionMap::GetImp().Clear();
}

void ActionMap::Free()
{
    // Historically deleted the ActionMapImp. There is no memory to
    // release now, so Free() clears all registered data instead,
    // preserving the "state resets after Free()" contract that existing
    // callers and tests depend on.
    ActionMap::GetImp().Clear();
}

int ActionMap::GetActionId( const std::string & name )
{
    return ActionMap::GetImp().GetActionId( name );
}

std::string ActionMap::GetActionName( int id )
{
    return ActionMap::GetImp().GetActionName( id );
}

void ActionMap::Register( const std::string & name )
{
    ActionMap::GetImp().Register( name );
}

void ActionMap::ReadFile( const std::string & fileName )
{
    ActionMap::GetImp().ReadFile( fileName );
}

EndNameSpace

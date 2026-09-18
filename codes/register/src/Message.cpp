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

#include "Message.h"
#include "TextFileParser.h"

BeginNameSpace( ONEFLOW )

void MessageMapImp::Register( const std::string & msgName )
{
    // Same defensive guard as ActionMapImp::Register: an empty name is
    // never a valid message and must not silently occupy an id slot.
    if ( msgName.empty() )
    {
        return;
    }

    if ( this->nameToId.find( msgName ) != this->nameToId.end() )
    {
        return; // first registration wins, unchanged from original behavior
    }

    int msgId = static_cast< int >( this->idToName.size() );
    this->nameToId[ msgName ] = msgId;
    this->idToName.push_back( msgName );
}

void MessageMapImp::Unregister( const std::string & msgName )
{
    // NOTE: same pre-existing asymmetry as ActionMapImp::Unregister -
    // this only removes the name->id mapping; the id->name slot in
    // idToName is left in place. Flagged, not fixed here, to keep this
    // step's scope limited to what MessageMap and ActionMap already had
    // in common. No caller currently exercises Unregister().
    this->nameToId.erase( msgName );
}

int MessageMapImp::GetMsgId( const std::string & msgName ) const
{
    auto iter = this->nameToId.find( msgName );
    if ( iter == this->nameToId.end() )
    {
        return -1;
    }
    return iter->second;
}

const std::string & MessageMapImp::GetMsgName( int msgId ) const
{
    static const std::string empty;
    if ( msgId < 0 || msgId >= static_cast< int >( this->idToName.size() ) )
    {
        return empty;
    }
    return this->idToName[ msgId ];
}

void MessageMapImp::ReadFile( const std::string & fileName )
{
    std::string separator = " =\r\n\t#$,;\"";

    TextFileParser textFileParser;
    textFileParser.OpenFile( fileName, std::ios_base::in );
    textFileParser.SetDefaultSeparator( separator );

    // FIX: was previously driven by ReachTheEndOfFile() + ReadNextNonEmptyLine(),
    // which had two bugs identical to the ones found and fixed in
    // ActionMapImp::ReadFile:
    //   1. eof() is only set AFTER a failed read, so the loop ran one
    //      extra time past the last line and registered a stray "".
    //   2. ReadNextNonEmptyLine() only skips blank lines, not comment
    //      lines ('#' or '//'), so a comment line's first token got
    //      registered as if it were real data.
    // Driving the loop by ReadNextMeaningfulLine()'s return value fixes
    // both: it skips blanks AND comments, and reliably reports "no more
    // content" at EOF.
    while ( textFileParser.ReadNextMeaningfulLine() )
    {
        std::string msgName = textFileParser.ReadNextWord();
        this->Register( msgName );
    }

    textFileParser.CloseFile();
}

void MessageMapImp::Clear()
{
    this->nameToId.clear();
    this->idToName.clear();
    ++ this->epoch_;
}

MessageMapImp & MessageMap::GetImp()
{
    // Meyer's singleton: same rationale as ActionMap::GetImp() - removes
    // the entire class of bugs (use-before-Init, dangling pointer after
    // Free, double-free) by construction, with no manual new/delete.
    static MessageMapImp imp;
    return imp;
}

void MessageMap::Init()
{
    MessageMap::GetImp().Clear();
}

void MessageMap::Free()
{
    MessageMap::GetImp().Clear();
}

int MessageMap::GetMsgId( const std::string & msgName )
{
    return MessageMap::GetImp().GetMsgId( msgName );
}

const std::string & MessageMap::GetMsgName( int msgId )
{
    return MessageMap::GetImp().GetMsgName( msgId );
}

int MessageMap::Epoch()
{
    return MessageMap::GetImp().Epoch();
}

void MessageMap::Register( const std::string & msgName )
{
    MessageMap::GetImp().Register( msgName );
}

bool MessageMap::Contains( const std::string & msgName )
{
    return MessageMap::GetMsgId( msgName ) >= 0;
}

void MessageMap::ReadFile( const std::string & fileName )
{
    MessageMap::GetImp().ReadFile( fileName );
}

EndNameSpace
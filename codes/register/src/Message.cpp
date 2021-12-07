/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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
#include "FileIO.h"

BeginNameSpace( ONEFLOW )
map< std::string, int > * MessageMap::nameMap = 0;
map< int, std::string > * MessageMap::idMap = 0;

MessageMap::MessageMap()
{
}

MessageMap::~MessageMap()
{
}

void MessageMap::Init()
{
    if ( MessageMap::nameMap ) return;
    MessageMap::nameMap = new std::map< std::string, int >();
    MessageMap::idMap = new std::map< int, std::string >();
}

void MessageMap::Free()
{
    delete MessageMap::nameMap;
    delete MessageMap::idMap;
    MessageMap::nameMap = 0;
    MessageMap::idMap = 0;
}

void MessageMap::Register( const std::string & msgName )
{
    std::map< std::string, int >::iterator iter = MessageMap::nameMap->find( msgName );
    if ( iter == MessageMap::nameMap->end() )
    {
        int msgId = MessageMap::nameMap->size();
        ( * MessageMap::nameMap )[ msgName ] = msgId;
        ( * MessageMap::idMap   )[ msgId   ] = msgName;
    }
}

void MessageMap::Unregister( const std::string & msgName )
{
    MessageMap::nameMap->erase( msgName );
}

int MessageMap::GetMsgId( const std::string & msgName )
{
    std::map< std::string, int >::iterator iter = MessageMap::nameMap->find( msgName );
    if ( iter == MessageMap::nameMap->end() )
    {
        return -1;
    }

    int actionIndex = iter->second;
    return actionIndex;
}

string MessageMap::GetMsgName( int msgId )
{
    std::map< int, std::string >::iterator iter = MessageMap::idMap->find( msgId );
    if ( iter == MessageMap::idMap->end() )
    {
        return "";
    }

    return iter->second;
}

void MessageMap::ReadFile( const std::string & fileName )
{
    std::string word;

    //\t is the tab key
    std::string separator = " =\r\n\t#$,;\"";

    FileIO ioFile;
    ioFile.OpenFile( fileName, std::ios_base::in );
    ioFile.SetDefaultSeparator( separator );

    while ( ! ioFile.ReachTheEndOfFile() )
    {
        ioFile.ReadNextNonEmptyLine();
        std::string msgName = ioFile.ReadNextWord();
        MessageMap::Register( msgName );
    }

    ioFile.CloseFile();
}

EndNameSpace

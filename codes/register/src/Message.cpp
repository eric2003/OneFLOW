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

#include "Message.h"
#include "FileIO.h"

BeginNameSpace( ONEFLOW )
map< string, int > * MessageMap::nameMap = 0;
map< int, string > * MessageMap::idMap = 0;

MessageMap::MessageMap()
{
}

MessageMap::~MessageMap()
{
}

void MessageMap::Init()
{
    if ( MessageMap::nameMap ) return;
    MessageMap::nameMap = new map< string, int >();
    MessageMap::idMap = new map< int, string >();
}

void MessageMap::Free()
{
    delete MessageMap::nameMap;
    delete MessageMap::idMap;
    MessageMap::nameMap = 0;
    MessageMap::idMap = 0;
}

void MessageMap::Register( const string & msgName )
{
    map< string, int >::iterator iter = MessageMap::nameMap->find( msgName );
    if ( iter == MessageMap::nameMap->end() )
    {
        int msgId = MessageMap::nameMap->size();
        ( * MessageMap::nameMap )[ msgName ] = msgId;
        ( * MessageMap::idMap   )[ msgId   ] = msgName;
    }
}

void MessageMap::Unregister( const string & msgName )
{
    MessageMap::nameMap->erase( msgName );
}

int MessageMap::GetMsgId( const string & msgName )
{
    map< string, int >::iterator iter = MessageMap::nameMap->find( msgName );
    if ( iter == MessageMap::nameMap->end() )
    {
        return -1;
    }

    int actionIndex = iter->second;
    return actionIndex;
}

string MessageMap::GetMsgName( int msgId )
{
    map< int, string >::iterator iter = MessageMap::idMap->find( msgId );
    if ( iter == MessageMap::idMap->end() )
    {
        return "";
    }

    return iter->second;
}

void MessageMap::ReadFile( const string & fileName )
{
    string word;

    //\tÎªtab¼ü
    string separator = " =\r\n\t#$,;\"";

    FileIO ioFile;
    ioFile.OpenFile( fileName, ios_base::in );
    ioFile.SetDefaultSeparator( separator );

    while ( ! ioFile.ReachTheEndOfFile() )
    {
        ioFile.ReadNextNonEmptyLine();
        string msgName = ioFile.ReadNextWord();
        MessageMap::Register( msgName );
    }

    ioFile.CloseFile();
}

EndNameSpace
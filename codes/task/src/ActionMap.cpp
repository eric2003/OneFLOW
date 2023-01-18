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
#include "ActionMap.h"
#include "DataBase.h"
#include "FileIO.h"

BeginNameSpace( ONEFLOW )

ActionMapImp * ActionMap::imp = 0;

ActionMap::ActionMap()
{
    ;
}

ActionMap::~ActionMap()
{
    ;
}

void ActionMap::Init()
{
    if ( ActionMap::imp ) return;
    ActionMap::imp = new ActionMapImp();
}

void ActionMap::Free()
{
    delete ActionMap::imp;
    ActionMap::imp = 0;
}

int ActionMap::GetActionId( const std::string & name )
{
    return ActionMap::imp->GetActionId( name );
}

std::string ActionMap::GetActionName( int id )
{
    return ActionMap::imp->GetActionName( id );
}

void ActionMap::ReadFile( const std::string & fileName )
{
    ActionMap::imp->ReadFile( fileName );
}

ActionMapImp::ActionMapImp()
{
    this->nameMap = new std::map< std::string, int >();
    this->idMap = new std::map< int, std::string >();
}

ActionMapImp::~ActionMapImp()
{
    delete this->nameMap;
    delete this->idMap;
}

void ActionMapImp::Register( const std::string & actionName )
{
    std::map< std::string, int >::iterator iter = this->nameMap->find( actionName );
    if ( iter == this->nameMap->end() )
    {
        int actionId = this->nameMap->size();
        ( * this->nameMap )[ actionName ] = actionId;
        ( * this->idMap   )[ actionId ] = actionName;
    }
}

void ActionMapImp::Unregister( const std::string & actionName )
{
    this->nameMap->erase( actionName );
}

int ActionMapImp::GetActionId( const std::string & actionName )
{
    std::map< std::string, int >::iterator iter = this->nameMap->find( actionName );
    if ( iter == this->nameMap->end() )
    {
        return -1;
    }

    int actionId = iter->second;
    return actionId;
}

std::string ActionMapImp::GetActionName( int actionIndex )
{
    std::map< int, std::string >::iterator iter = this->idMap->find( actionIndex );
    if ( iter == this->idMap->end() )
    {
        return "";
    }

    std::string & actionName = iter->second;
    return actionName;
}

void ActionMapImp::ReadFile( const std::string & fileName )
{
    //\t is the tab key
    std::string separator = " =\r\n\t#$,;\"";

    FileIO ioFile;
    ioFile.OpenFile( fileName, std::ios_base::in );
    ioFile.SetDefaultSeparator( separator );

    while ( ! ioFile.ReachTheEndOfFile() )
    {
        ioFile.ReadNextNonEmptyLine();
        std::string actionName = ioFile.ReadNextWord();
        this->Register( actionName );
    }

    ioFile.CloseFile();
}



EndNameSpace

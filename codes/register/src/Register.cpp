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

#include "Register.h"
#include "RegisterUtil.h"
#include "Category.h"
#include "SolverInfo.h"
#include "HXClone.h"
#include "FileIO.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )
HXRegister::HXRegister()
{
    ;
}

HXRegister::~HXRegister()
{
    ;
}

void HXRegister::Register( const string & cmdName, const string & className )
{
    map< string, HXClone * >::iterator iter = this->data.find( cmdName );
    if ( iter == this->data.end() )
    {
        HXClone * cloneClass = HXClone::SafeClone( className );
        this->data[ cmdName ] = cloneClass;
    }
}

HXClone * HXRegister::GetClass( const string & cmdName )
{
    map< string, HXClone * >::iterator iter = this->data.find( cmdName );

    if ( iter != this->data.end() )
    {
        return iter->second;
    }
    return 0;
}

void HXRegister::FreeAll()
{
    for ( auto iter = this->data.begin(); iter != this->data.end(); ++ iter )
    {
        delete iter->second;
    }
    this->data.clear();
}

MRegister::MRegister()
{
}

MRegister::~MRegister()
{
    int nRegisters = this->data.size();
    for ( int iRegister = 0; iRegister < nRegisters; ++ iRegister )
    {
        delete this->data[ iRegister ];
    }
}

HXRegister * MRegister::GetRegister( int index )
{
    return this->data[ index ];
}

HXRegister * MRegister::GetRegister()
{
    int index = 0;
    return this->data[ index ];
}

void MRegister::AllocateData()
{
    int nRegisters = this->fileNames.size();

    if ( this->data.size() == nRegisters ) return;

    for ( int iRegister = 0; iRegister < nRegisters; ++ iRegister )
    {
        this->data.push_back( new HXRegister() );
    }
}

void MRegister::SetSolverFileNames( StringField & fileNames )
{
    this->fileNames = fileNames;
}

void MRegister::RegisterAll()
{
    this->AllocateData();

    for ( int iRegister = 0; iRegister < this->data.size(); ++ iRegister )
    {
        HXRegister * fRegister = this->data[ iRegister ];
        string & fileName = fileNames[ iRegister ];
        this->Register( fileName, fRegister );
    }
}

void MRegister::Register( const string & fileName, HXRegister * fRegister )
{
    //\tÎªtab¼ü
    string separator  = " =\r\n\t#$,;\"()";

    FileIO ioFile;
    ioFile.OpenFile( fileName, ios_base::in );
    ioFile.SetDefaultSeparator( separator );

    while ( ! ioFile.ReachTheEndOfFile()  )
    {
        bool flag = ioFile.ReadNextNonEmptyLine();
        if ( ! flag ) break;
        string actionName = ioFile.ReadNextWord();
        string className  = ioFile.ReadNextWord();
        //cout << "actionName = " << actionName << " className= " << className << endl;
        if ( className == "CFillWallStruct" )
        {
            int kkk = 1;
        }
        fRegister->Register( actionName, className );

        HXClone * cloneClass = fRegister->GetClass( actionName );
        int nParameters = ioFile.ReadNextDigit< int >();
        for ( int iParameter = 0; iParameter < nParameters; ++ iParameter )
        {
            cloneClass->data.push_back( ioFile.ReadNextWord() );
        }
    }

    ioFile.CloseFile();
}

map< int, MRegister * > * RegisterFactory::data = 0;

RegisterFactory::RegisterFactory()
{
}

RegisterFactory::~RegisterFactory()
{
}

void RegisterFactory::Init()
{
    if ( ! RegisterFactory::data )
    {
        RegisterFactory::data = new map< int, MRegister * >();
    }
}

void RegisterFactory::AddMRegister( int registerId )
{
    map< int, MRegister * >::iterator iter;
    RegisterFactory::Init();
    iter = RegisterFactory::data->find( registerId );
    if ( iter == RegisterFactory::data->end() )
    {
        MRegister * mRegister = new MRegister();
        ( * RegisterFactory::data )[ registerId ] = mRegister;
    }
}

MRegister * RegisterFactory::GetMRegister( int registerId )
{
    map< int, MRegister * >::iterator iter;
    iter = RegisterFactory::data->find( registerId );
    return iter->second;
}

void RegisterFactory::FreeMRegister()
{
    if ( ! RegisterFactory::data ) return;
    map< int, MRegister * >::iterator iter;
    for ( iter = RegisterFactory::data->begin(); iter != RegisterFactory::data->end(); ++ iter )
    {
        delete iter->second;
    }

    RegisterFactory::data->clear();

    delete RegisterFactory::data;
    RegisterFactory::data = 0;
}

HXRegister * RegisterFactory::GetRegister( int mRegisterId, int registerId )
{
    MRegister * mRegister = RegisterFactory::GetMRegister( mRegisterId );
    HXRegister * hRegister = mRegister->GetRegister( registerId );
    return hRegister;
}

EndNameSpace
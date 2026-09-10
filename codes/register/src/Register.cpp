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

//#include "Register.h"
//#include "RegisterUtils.h"
//#include "Category.h"
//#include "SolverInfo.h"
//#include "HXClone.h"
//#include "TextFileParser.h"
//#include <iostream>
//
//
//BeginNameSpace( ONEFLOW )
//HXRegister::HXRegister()
//{
//    ;
//}
//
//HXRegister::~HXRegister()
//{
//    ;
//}
//
//void HXRegister::Register( const std::string & cmdName, const std::string & className )
//{
//    std::map< std::string, HXClone * >::iterator iter = this->data.find( cmdName );
//    if ( iter == this->data.end() )
//    {
//        HXClone * cloneClass = HXClone::SafeClone( className );
//        this->data[ cmdName ] = cloneClass;
//    }
//}
//
//HXClone * HXRegister::GetClass( const std::string & cmdName )
//{
//    std::map< std::string, HXClone * >::iterator iter = this->data.find( cmdName );
//
//    if ( iter != this->data.end() )
//    {
//        return iter->second;
//    }
//    return 0;
//}
//
//void HXRegister::FreeAll()
//{
//    for ( auto iter = this->data.begin(); iter != this->data.end(); ++ iter )
//    {
//        delete iter->second;
//    }
//    this->data.clear();
//}
//
//MRegister::MRegister()
//{
//}
//
//MRegister::~MRegister()
//{
//    int nRegisters = this->data.size();
//    for ( int iRegister = 0; iRegister < nRegisters; ++ iRegister )
//    {
//        delete this->data[ iRegister ];
//    }
//}
//
//HXRegister * MRegister::GetRegister( int index )
//{
//    return this->data[ index ];
//}
//
//HXRegister * MRegister::GetRegister()
//{
//    int index = 0;
//    return this->data[ index ];
//}
//
//void MRegister::AllocateData()
//{
//    int nRegisters = this->fileNames.size();
//
//    if ( this->data.size() == nRegisters ) return;
//
//    for ( int iRegister = 0; iRegister < nRegisters; ++ iRegister )
//    {
//        this->data.push_back( new HXRegister() );
//    }
//}
//
//void MRegister::SetSolverFileNames( StringField & fileNames )
//{
//    this->fileNames = fileNames;
//}
//
//void MRegister::RegisterAll()
//{
//    this->AllocateData();
//
//    for ( int iRegister = 0; iRegister < this->data.size(); ++ iRegister )
//    {
//        HXRegister * fRegister = this->data[ iRegister ];
//        std::string & fileName = fileNames[ iRegister ];
//        this->Register( fileName, fRegister );
//    }
//}
//
//void MRegister::Register( const std::string & fileName, HXRegister * fRegister )
//{
//    //\t is the tab key
//    std::string separator  = " =\r\n\t#$,;\"()";
//
//    TextFileParser textFileParser;
//    textFileParser.OpenFile( fileName, std::ios_base::in );
//    textFileParser.SetDefaultSeparator( separator );
//
//    while ( ! textFileParser.ReachTheEndOfFile()  )
//    {
//        bool flag = textFileParser.ReadNextNonEmptyLine();
//        if ( ! flag ) break;
//        std::string actionName = textFileParser.ReadNextWord();
//        std::string className  = textFileParser.ReadNextWord();
//        //std::cout << "actionName = " << actionName << " className= " << className << std::endl;
//
//        fRegister->Register( actionName, className );
//
//        HXClone * cloneClass = fRegister->GetClass( actionName );
//        int nParameters = textFileParser.ReadNextDigit< int >();
//        for ( int iParameter = 0; iParameter < nParameters; ++ iParameter )
//        {
//            cloneClass->data.push_back( textFileParser.ReadNextWord() );
//        }
//    }
//
//    textFileParser.CloseFile();
//}
//
//std::map< int, MRegister * > * RegisterFactory::data = 0;
//
//RegisterFactory::RegisterFactory()
//{
//}
//
//RegisterFactory::~RegisterFactory()
//{
//}
//
//void RegisterFactory::Init()
//{
//    if ( ! RegisterFactory::data )
//    {
//        RegisterFactory::data = new std::map< int, MRegister * >();
//    }
//}
//
//void RegisterFactory::AddMRegister( int registerId )
//{
//    std::map< int, MRegister * >::iterator iter;
//    RegisterFactory::Init();
//    iter = RegisterFactory::data->find( registerId );
//    if ( iter == RegisterFactory::data->end() )
//    {
//        MRegister * mRegister = new MRegister();
//        ( * RegisterFactory::data )[ registerId ] = mRegister;
//    }
//}
//
//MRegister * RegisterFactory::GetMRegister( int registerId )
//{
//    std::map< int, MRegister * >::iterator iter;
//    iter = RegisterFactory::data->find( registerId );
//    return iter->second;
//}
//
//void RegisterFactory::FreeMRegister()
//{
//    if ( ! RegisterFactory::data ) return;
//    std::map< int, MRegister * >::iterator iter;
//    for ( iter = RegisterFactory::data->begin(); iter != RegisterFactory::data->end(); ++ iter )
//    {
//        delete iter->second;
//    }
//
//    RegisterFactory::data->clear();
//
//    delete RegisterFactory::data;
//    RegisterFactory::data = 0;
//}
//
//HXRegister * RegisterFactory::GetRegister( int mRegisterId, int registerId )
//{
//    MRegister * mRegister = RegisterFactory::GetMRegister( mRegisterId );
//    HXRegister * hRegister = mRegister->GetRegister( registerId );
//    return hRegister;
//}
//
//EndNameSpace


#include "Register.h"
#include "RegisterUtils.h"
#include "Category.h"
#include "SolverInfo.h"
#include "HXClone.h"
#include "TextFileParser.h"
#include <iostream>

BeginNameSpace( ONEFLOW )

void HXRegister::Register( const std::string & cmdName, const std::string & className )
{
    if ( this->data.find( cmdName ) != this->data.end() )
    {
        return; // already registered, unchanged from original behavior
    }

    HXClone * cloneClass = HXClone::SafeClone( className );
    this->data[ cmdName ] = std::unique_ptr< HXClone >( cloneClass );
}

HXClone * HXRegister::GetClass( const std::string & cmdName )
{
    auto iter = this->data.find( cmdName );
    if ( iter != this->data.end() )
    {
        return iter->second.get();
    }
    return nullptr;
}

void HXRegister::FreeAll()
{
    // Clearing the map runs each unique_ptr's destructor, releasing
    // every owned HXClone. Kept as a named method for source
    // compatibility with any existing explicit FreeAll() call sites.
    this->data.clear();
}

HXRegister * MRegister::GetRegister( int index )
{
    // FIX: previously `return this->data[index];` ¡ª out-of-bounds access
    // via vector::operator[] is undefined behavior. Bounds-check and
    // return nullptr instead; callers must handle the nullptr case.
    if ( index < 0 || index >= static_cast< int >( this->data.size() ) )
    {
        return nullptr;
    }
    return this->data[ index ].get();
}

HXRegister * MRegister::GetRegister()
{
    return this->GetRegister( 0 );
}

void MRegister::AllocateData()
{
    int nRegisters = static_cast< int >( this->fileNames.size() );

    if ( static_cast< int >( this->data.size() ) == nRegisters ) return;

    for ( int iRegister = 0; iRegister < nRegisters; ++ iRegister )
    {
        this->data.push_back( std::make_unique< HXRegister >() );
    }
}

void MRegister::SetSolverFileNames( StringField & fileNames )
{
    this->fileNames = fileNames;
}

void MRegister::RegisterAll()
{
    this->AllocateData();

    for ( int iRegister = 0; iRegister < static_cast< int >( this->data.size() ); ++ iRegister )
    {
        HXRegister * fRegister = this->data[ iRegister ].get();
        std::string & fileName = fileNames[ iRegister ];
        this->Register( fileName, fRegister );
    }
}

void MRegister::Register( const std::string & fileName, HXRegister * fRegister )
{
    std::string separator = " =\r\n\t#$,;\"()";

    TextFileParser textFileParser;
    textFileParser.OpenFile( fileName, std::ios_base::in );
    textFileParser.SetDefaultSeparator( separator );

    // NOTE: left as ReachTheEndOfFile()+ReadNextNonEmptyLine() for now -
    // this loop also reads a parameter count and N following words per
    // entry, which is a different and more complex shape than the
    // "one name per line" pattern fixed elsewhere. Changing its
    // end-of-file handling deserves its own dedicated look rather than
    // a drive-by change bundled into this ownership refactor.
    while ( ! textFileParser.ReachTheEndOfFile() )
    {
        bool flag = textFileParser.ReadNextNonEmptyLine();
        if ( ! flag ) break;
        std::string actionName = textFileParser.ReadNextWord();
        std::string className  = textFileParser.ReadNextWord();

        fRegister->Register( actionName, className );

        HXClone * cloneClass = fRegister->GetClass( actionName );
        int nParameters = textFileParser.ReadNextDigit< int >();
        for ( int iParameter = 0; iParameter < nParameters; ++ iParameter )
        {
            cloneClass->data.push_back( textFileParser.ReadNextWord() );
        }
    }

    textFileParser.CloseFile();
}

std::map< int, std::unique_ptr< MRegister > > & RegisterFactory::GetData()
{
    // Meyer's singleton, same rationale as ActionMap/MessageMap: no
    // manual new/delete, so use-before-Init and double-free are
    // structurally impossible.
    static std::map< int, std::unique_ptr< MRegister > > data;
    return data;
}

void RegisterFactory::Init()
{
    // Historically allocated the map itself. GetData() now self-
    // initializes, so Init() is kept only for source compatibility and
    // is a safe no-op to call any number of times.
}

void RegisterFactory::AddMRegister( int registerId )
{
    auto & data = RegisterFactory::GetData();
    if ( data.find( registerId ) == data.end() )
    {
        data[ registerId ] = std::make_unique< MRegister >();
    }
}

MRegister * RegisterFactory::GetMRegister( int registerId )
{
    // FIX: previously `return iter->second;` without checking for
    // end() - dereferencing end() on an unregistered registerId was
    // undefined behavior. Return nullptr instead; callers must check.
    auto & data = RegisterFactory::GetData();
    auto iter = data.find( registerId );
    if ( iter == data.end() )
    {
        return nullptr;
    }
    return iter->second.get();
}

void RegisterFactory::FreeMRegister()
{
    // Clearing the map runs each unique_ptr<MRegister>'s destructor,
    // which in turn runs each unique_ptr<HXRegister>'s destructor, which
    // in turn releases every owned HXClone. The entire ownership chain
    // is now cleaned up automatically - no more manual delete loop, and
    // no more silent leak of the innermost HXClone objects.
    RegisterFactory::GetData().clear();
}

HXRegister * RegisterFactory::GetRegister( int mRegisterId, int registerId )
{
    // FIX: previously dereferenced GetMRegister's result unconditionally.
    // Now short-circuits to nullptr if either lookup fails.
    MRegister * mRegister = RegisterFactory::GetMRegister( mRegisterId );
    if ( ! mRegister )
    {
        return nullptr;
    }
    return mRegister->GetRegister( registerId );
}

EndNameSpace

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

#include "DataPage.h"
#include "BasicParallel.h"
#include "Parallel.h"
#include "Stop.h"

#ifndef _WINDOWS
   #include <string.h>
#endif
#include <fstream>


BeginNameSpace( ONEFLOW )
DataPage::DataPage()
{
    this->dataMemory = new CharMemory();
    this->currPos = 0;
}

DataPage::~DataPage()
{
    delete this->dataMemory;
}

void DataPage::MoveToPosition( HXSize_t position )
{
    if ( ( 0 < position ) && ( position < GetSize() ) )
    {
        this->currPos = position;
    }
    else if ( 0 == position )
    {
        this->currPos = position;
    }
    else
    {
        Stop( "Out of Range: position \n" );
    }
}

void DataPage::MoveForwardPosition( HXSize_t dataSize )
{
    this->currPos += dataSize;
}

HXSize_t DataPage::GetSize()
{
    return dataMemory->size();
}

char * DataPage::GetBeginDataPointer()
{
    if ( this->GetSize() == 0 )
    {
        return 0;
    }
    return &( ( * dataMemory )[ 0 ] );
}

char * DataPage::GetCurrentDataPointer()
{
    if ( this->GetSize() == 0 )
    {
        return 0;
    }
    return &( ( * dataMemory )[ currPos ] );
}

void DataPage::ToString( std::string & str )
{
    if ( this->GetSize() )
    {
        str.append( this->GetBeginDataPointer(), this->GetSize() );
    }
}

void DataPage::Write( void * data, HXSize_t dataSize )
{
    if ( dataSize <= 0 ) return;

    memcpy( this->GetCurrentDataPointer(), data, dataSize );
    this->MoveForwardPosition( dataSize );
}

void DataPage::Read( void * data, HXSize_t dataSize )
{
    if ( dataSize <= 0 ) return;
    memcpy( data, this->GetCurrentDataPointer(), dataSize );
    this->MoveForwardPosition( dataSize );
}

void DataPage::Write( void * data, HXSize_t dataSize, HXSize_t position )
{
    this->MoveToPosition( position );
    this->Write( data, dataSize );
}

void DataPage::Read( void * data, HXSize_t dataSize, HXSize_t position )
{
    this->MoveToPosition( position );
    this->Read( data, dataSize );
}

void DataPage::ReSize( HXSize_t newSize )
{
    this->dataMemory->resize( newSize );
}

void DataPage::Send( int pId, int tag )
{
    HXSize_t nLength = this->GetSize();

    if ( nLength <= 0 ) return;
    ONEFLOW::HXSend( this->GetBeginDataPointer(), nLength, PL_CHAR, pId, tag );
}

void DataPage::Recv( int pId, int tag )
{
    HXSize_t nLength = this->GetSize();

    if ( nLength <= 0 ) return;

    ONEFLOW::HXRecv( this->GetBeginDataPointer(), nLength, PL_CHAR, pId, tag );
}

void DataPage::Bcast( int rootid )
{
    HXSize_t nLength = this->GetSize();

    if ( nLength <= 0 ) return;
    HXBcast( this->GetBeginDataPointer(), nLength, rootid );
}

void DataPage::ReadFile( std::fstream & file )
{
    HXSize_t nLength = this->GetSize();

    if ( nLength <= 0 ) return;

    char * data = this->GetBeginDataPointer();

    file.read( data, nLength );
}

void DataPage::WriteFile( std::fstream & file )
{
    HXSize_t nLength = this->GetSize();

    if ( nLength <= 0 )
    {
        return;
    }

    char * data = this->GetBeginDataPointer();

    file.write( data, nLength );
}

EndNameSpace

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

#include "DataPage.h"
#include "BasicParallel.h"
#include "Parallel.h"
#include "Fatal.h"
#include <cstring>

//#ifndef _WINDOWS
//   #include <string.h>
//#endif
#include <fstream>


BeginNameSpace( ONEFLOW )
DataPage::DataPage()
{
    this->currPos = 0;
}

DataPage::~DataPage()
{
}

void DataPage::MoveToPosition( HXSize_t position )
{
    // position == GetSize() is a valid "end/append" position.
    if ( position <= size() )
    {
        this->currPos = position;
    }
    else
    {
        Fatal( "Out of Range: position \n" );
    }
}

void DataPage::Advance( HXOffset_t offset )
{
    this->currPos += offset;
}

HXSize_t DataPage::size() const
{
    return dataMemory.size();
}

char * DataPage::data()
{
    // Use data() instead of operator[] to avoid UB.
    return dataMemory.data();
}

const char * DataPage::data() const
{
    // Provide a const overload for read-only access.
    return dataMemory.data();
}

char * DataPage::CurrentPtr()
{
    // currPos may legitimately equal GetSize() (the "append/end" position).
    // vector::data() + size() is well-defined as long as it is never dereferenced,
    // unlike operator[](size()) which is UB even just to take its address.
    //return dataMemory.data() + currPos;
    return this->PtrAt( currPos );
}

char * DataPage::PtrAt( int offset )
{
    // Same reasoning as above; offset == size() must be safe to compute.
    return dataMemory.data() + offset;
}

void DataPage::ToString( std::string & str )
{
    if ( this->size() )
    {
        str.append( this->data(), this->size() );
    }
}

void DataPage::Write( const void * data, HXSize_t dataSize )
{
    if ( dataSize == 0 || data == nullptr )
    {
        return;
    }

    // Check for buffer overflow before writing.
    if ( this->currPos + dataSize > this->size() )
    {
        throw std::out_of_range("DataPage::Write - Buffer overflow");
    }

    std::memcpy( this->CurrentPtr(), data, dataSize );
    this->Advance( dataSize );
}

void DataPage::Read( void * data, HXSize_t dataSize )
{
    if ( dataSize == 0 || data == nullptr )
    {
        return;
    }

    // Prevent potential integer overflow and check bound.
    if ( dataSize > this->size() - this->currPos )
    {
        throw std::out_of_range("DataPage::Read - Attempted to read past end of buffer");
    }

    std::memcpy( data, this->CurrentPtr(), dataSize );
    this->Advance( dataSize );
}

void DataPage::Write( const void * data, HXSize_t position, HXSize_t dataSize )
{
    if ( dataSize == 0 || data == nullptr )
    {
        return;
    }

    // Check bounds for the specific window without changing internal state.
    if ( position > this->size() || dataSize > this->size() - position )
    {
        throw std::out_of_range("DataPage::Write - Buffer overflow at specified position");
    }

    std::memcpy( this->data() + position, data, dataSize );
}

void DataPage::Read( void * data, HXSize_t position, HXSize_t dataSize ) const
{
    if ( dataSize == 0 || data == nullptr )
    {
        return;
    }

    // Check bounds for the specific window without changing internal state.
    if ( position > this->size() || dataSize > this->size() - position )
    {
        throw std::out_of_range("DataPage::Read - Attempted to read past end of buffer");
    }

    std::memcpy( data, this->data() + position, dataSize );
}

void DataPage::ReSize( HXSize_t newSize )
{
    this->dataMemory.resize( newSize );
}

void DataPage::Send( int pId, int tag )
{
    HXSize_t nLength = this->size();

    if ( nLength <= 0 ) return;
    ONEFLOW::HXSend( this->data(), nLength, PL_CHAR, pId, tag );
}

void DataPage::Recv( int pId, int tag )
{
    HXSize_t nLength = this->size();

    if ( nLength <= 0 ) return;

    ONEFLOW::HXRecv( this->data(), nLength, PL_CHAR, pId, tag );
}

void DataPage::Bcast( int rootid )
{
    HXSize_t nLength = this->size();

    if ( nLength <= 0 ) return;
    HXBcast( this->data(), nLength, rootid );
}

void DataPage::ReadFile( std::fstream & file )
{
    HXSize_t nLength = this->size();

    if ( nLength <= 0 ) return;

    char * data = this->data();

    file.read( data, nLength );
}

void DataPage::WriteFile( std::fstream & file )
{
    HXSize_t nLength = this->size();

    if ( nLength <= 0 )
    {
        return;
    }

    char * data = this->data();

    file.write( data, nLength );
}

EndNameSpace

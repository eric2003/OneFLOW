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


#include "DataBook.h"
#include "DataBaseIO.h"
#include "BasicParallel.h"
#include "Parallel.h"
#include <sstream>
using namespace std;

BeginNameSpace( ONEFLOW )

DataBook::DataBook()
{
    //maxUnitSize = 1024000;
    maxUnitSize = 1024000000;
    dataBook = new vector< DataPage * >;
    dataBook->push_back( new DataPage() );
    this->currPos = 0;
    this->currPageId = 0;

    this->MoveToBegin();
}

DataBook::~DataBook()
{
    for ( UInt i = 0; i < dataBook->size(); ++ i )
    {
        delete ( * dataBook )[ i ];
    }
    delete dataBook;
}

UInt DataBook::GetNPage()
{
    return dataBook->size();
}

DataPage * DataBook::GetCurrentPage()
{
    currPageId = this->currPos / maxUnitSize;
    return ( * dataBook )[ currPageId ];
}

DataPage * DataBook::GetPage( UInt iPage )
{
    return ( * dataBook )[ iPage ];
}

char * MovePointer( void * data, LLong dataSize )
{
    return reinterpret_cast< char * >( data ) + dataSize;
}

void DataBook::MoveForwardPosition( LLong dataSize )
{
    this->currPos += dataSize;
}

void DataBook::Read( void * data, LLong dataSize )
{
    if ( dataSize <= 0 ) return;

    LLong remainingSize = this->GetRemainingSizeOfCurrentPage();

    if ( remainingSize >= dataSize )
    {
        this->GetCurrentPage()->Read( data, dataSize );
        this->MoveForwardPosition( dataSize );
    }
    else
    {
        this->GetCurrentPage()->Read( data, remainingSize );
        this->MoveForwardPosition( remainingSize );

        void * newData = ONEFLOW::MovePointer( data, remainingSize );
        LLong newSize = dataSize - remainingSize;

        this->Read( newData, newSize );
    }
}

void DataBook::Write( void * data, LLong dataSize )
{
    if ( dataSize <= 0 ) return;

    this->SecureRelativeSpace( dataSize );

    LLong remainingSize = this->GetRemainingSizeOfCurrentPage();

    if ( remainingSize >= dataSize )
    {
        this->GetCurrentPage()->Write( data, dataSize );
        this->MoveForwardPosition( dataSize );
    }
    else
    {
        this->GetCurrentPage()->Write( data, remainingSize );

        this->MoveForwardPosition( remainingSize );

        void * newData = ONEFLOW::MovePointer( data, remainingSize );

        LLong  newSize = dataSize - remainingSize;

        this->Write( newData, newSize );
    }
}

void DataBook::ReadString( string & cs )
{
    UInt nLength = 0;
    this->Read( & nLength, sizeof( UInt ) );

    char * data = new char[ nLength + 1 ];

    this->Read( data, nLength + 1 );

    cs = data;

    delete[] data;
}

void DataBook::WriteString( string & cs )
{
    UInt nLength = cs.length();

    this->Write( & nLength, sizeof( UInt ) );

    char * data = new char[ nLength + 1 ];

    cs.copy( data, nLength );
    data[ nLength ] = '\0';

    this->Write( data, nLength + 1 );

    delete[] data;
}

void DataBook::AppendString( string & cs )
{
    this->MoveToEnd();
    this->WriteString( cs );
}

void DataBook::Write( ostringstream * oss )
{
    string str = oss->str();
    UInt stringSize = str.size();
    this->Write( const_cast< char * >( str.c_str() ), stringSize * sizeof( char ) );
}

LLong DataBook::GetSize()
{
    LLong sum = 0;
    for ( int iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        sum += this->GetPage( iPage )->GetSize();
    }
    return sum;
}

void DataBook::ReSize( LLong nLength )
{
    if ( nLength <= 0 )
    {
        if ( nLength == 0 )
        {
            for ( UInt iPage = 0; iPage < this->GetNPage(); ++ iPage )
            {
                this->GetPage( iPage )->ReSize( 0 );
            }
        }
        return;
    }

    //23 divided by 3 is 7, remainder 2.
    UInt nPage = nLength / maxUnitSize;
    LLong remainder = nLength % maxUnitSize;

    UInt additionalPage = 0;
    if ( remainder )
    {
        additionalPage = 1;
    }

    UInt newNPage = nPage + additionalPage;
    this->ResizeNPage( newNPage );

    for ( UInt iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        LLong needSize = maxUnitSize;
        if ( iPage == nPage )
        {
            needSize = remainder;
        }
        this->GetPage( iPage )->ReSize( needSize );
    }
}

void DataBook::ResizeNPage( UInt newNPage )
{
    UInt oldNPage = this->GetNPage();

    if ( newNPage <= oldNPage )
    {
        this->Erase( newNPage, oldNPage );
        dataBook->resize( newNPage );
    }
    else
    {
        UInt iPageStart = oldNPage;
        UInt iPageEnd = newNPage;

        for ( UInt iPage = iPageStart; iPage != iPageEnd; ++ iPage )
        {
            dataBook->push_back( new DataPage() );
        }
    }
}

void DataBook::SecureRelativeSpace( LLong dataSize )
{
    LLong needSize = this->currPos + dataSize;

    this->SecureAbsoluteSpace( needSize );
}

void DataBook::SecureAbsoluteSpace( LLong needSize )
{
    //If there is enough space, there is no need to allocate
    //This can cause some string problems, and if not ReSize, there may be superfluous characters
    //in the memory that are not cleared
    //if ( needSize <= GetSize() ) return;

    this->ReSize( needSize );
}

void DataBook::MoveToBegin()
{
    this->currPos = 0;
    for ( int iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        this->GetPage( iPage )->MoveToBegin();
    }
}

void DataBook::MoveToEnd()
{
    this->currPos = this->GetSize();
    for ( UInt iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        this->GetPage( iPage )->MoveToEnd();
    }
}

LLong DataBook::GetRemainingSizeOfCurrentPage()
{
    LLong remainder = this->currPos % maxUnitSize;
    return maxUnitSize - remainder;
}

void DataBook::ReadFile( fstream & file )
{
    //Read the contents of file into DataBook
    //And for DataBook, the process is counter, equivalent to writing

    LLong nLength = 0;
    ONEFLOW::HXRead( & file, nLength );

    if ( nLength <= 0 ) return;

    this->SecureAbsoluteSpace( nLength );

    for ( streamsize iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        this->GetPage( iPage )->ReadFile( file );
    }
}

void DataBook::WriteFile( fstream & file )
{
    LLong nLength = this->GetSize();

    //Whether or not nLength is less than zero, you need to write the file
    ONEFLOW::HXWrite( & file, nLength );
    if ( nLength <= 0 )
    {
        return;
    }

    for ( UInt iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        this->GetPage( iPage )->WriteFile( file );
    }
}

void DataBook::ToString( string & str )
{
    for ( UInt iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        this->GetPage( iPage )->ToString( str );
    }
}

void DataBook::Append( void * data, LLong dataSize )
{
    this->MoveToEnd();
    LLong needSize = this->GetSize() + dataSize;
    this->SecureAbsoluteSpace( needSize );

    this->GetCurrentPage()->Write( data, dataSize );
}

void DataBook::Destroy( DataPage * dataPage )
{
    delete dataPage;
}

void DataBook::Erase( UInt startPage, UInt endPage )
{
    for ( UInt iPage = startPage; iPage != endPage; ++ iPage )
    {
        this->Destroy( GetPage( iPage ) );
    }
}

void DataBook::Send( int pid, int tag )
{
    LLong nLength = this->GetSize();

    ONEFLOW::HXSend( & nLength, 1, PL_LONG_LONG_INT, pid, tag );

    //It is necessary to judge the zero of data length
    if ( nLength <= 0 ) return;

    for ( UInt iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        this->GetPage( iPage )->Send( pid, tag );
    }
}

void DataBook::Recv( int pid, int tag )
{
    LLong nLength = 0;

    ONEFLOW::HXRecv( & nLength, 1, PL_LONG_LONG_INT, pid, tag );

    if ( nLength <= 0 )
    {
        return;
    }

    this->SecureAbsoluteSpace( nLength );

    for ( UInt iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        this->GetPage( iPage )->Recv( pid, tag );
    }
}

void DataBook::SendRecv( int sendpid, int recvpid, int tag )
{
    if ( sendpid == recvpid ) return;

    if ( Parallel::pid == sendpid )
    {
        this->Send( recvpid, tag );
    }
    else if ( Parallel::pid == recvpid )
    {
        this->Recv( sendpid, tag );
    }
}

void DataBook::Bcast( int rootid )
{
    LLong nLength = this->GetSize();

    HXBcast( & nLength, 1, rootid );

    if ( nLength <= 0 )
    {
        return;
    }

    if ( Parallel::pid != rootid )
    {
        this->SecureAbsoluteSpace( nLength );
    }

    for ( UInt iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        this->GetPage( iPage )->Bcast( rootid );
    }
}

void ToDataBook( DataBook * dataBook, ostringstream & oss )
{
    if ( ! dataBook ) return;

    dataBook->MoveToBegin();
    dataBook->ReSize( 0 );
    dataBook->Write( & oss );
}

EndNameSpace

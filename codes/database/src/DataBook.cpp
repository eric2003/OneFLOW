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


#include "DataBook.h"
#include "DataBaseIO.h"
#include "BasicParallel.h"
#include "Parallel.h"
#include <sstream>
#include "LogFile.h"

BeginNameSpace( ONEFLOW )

DataBook::DataBook()
{
    //maxUnitSize = 1024000;
    maxUnitSize = 1024000000;
    dataBook = new std::vector< DataPage * >;
    dataBook->push_back( new DataPage() );
    this->currPos = 0;
    this->currPageId = 0;

    this->MoveToBegin();
}

DataBook::~DataBook()
{
    for ( HXSize_t i = 0; i < dataBook->size(); ++ i )
    {
        delete ( * dataBook )[ i ];
    }
    delete dataBook;
}

HXSize_t DataBook::GetNPage()
{
    return dataBook->size();
}

DataPage * DataBook::GetCurrentPage()
{
    currPageId = this->currPos / maxUnitSize;
    return ( * dataBook )[ currPageId ];
}

DataPage * DataBook::GetPage( HXSize_t iPage )
{
    return ( * dataBook )[ iPage ];
}

char * MovePointer( void * data, HXLongLong_t dataSize )
{
    return reinterpret_cast< char * >( data ) + dataSize;
}

void DataBook::MoveForwardPosition( HXLongLong_t dataSize )
{
    this->currPos += dataSize;
}

void DataBook::Read( void * data, HXLongLong_t dataSize )
{
    if ( dataSize <= 0 ) return;

    HXLongLong_t remainingSize = this->GetRemainingSizeOfCurrentPage();

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
        HXLongLong_t newSize = dataSize - remainingSize;

        this->Read( newData, newSize );
    }
}

void DataBook::Write( void * data, HXLongLong_t dataSize )
{
    if ( dataSize <= 0 ) return;

    this->SecureRelativeSpace( dataSize );

    HXLongLong_t remainingSize = this->GetRemainingSizeOfCurrentPage();

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

        HXLongLong_t  newSize = dataSize - remainingSize;

        this->Write( newData, newSize );
    }
}

void DataBook::ReadString( std::string & cs )
{
    HXSize_t nLength = 0;
    this->Read( & nLength, sizeof( HXSize_t ) );

    char * data = new char[ nLength + 1 ];

    this->Read( data, nLength + 1 );

    cs = data;

    delete[] data;
}

void DataBook::WriteString( std::string & cs )
{
    HXSize_t nLength = cs.length();

    this->Write( & nLength, sizeof( HXSize_t ) );

    char * data = new char[ nLength + 1 ];

    cs.copy( data, nLength );
    data[ nLength ] = '\0';

    this->Write( data, nLength + 1 );

    delete[] data;
}

void DataBook::AppendString( std::string & cs )
{
    this->MoveToEnd();
    this->WriteString( cs );
}

void DataBook::Write( std::ostringstream * oss )
{
    std::string str = oss->str();
    HXSize_t stringSize = str.size();
    this->Write( const_cast< char * >( str.c_str() ), stringSize * sizeof( char ) );
}

HXLongLong_t DataBook::GetSize()
{
    HXLongLong_t sum = 0;
    for ( int iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        sum += this->GetPage( iPage )->GetSize();
    }
    return sum;
}

void DataBook::ReSize( HXLongLong_t nLength )
{
    if ( nLength <= 0 )
    {
        if ( nLength == 0 )
        {
            for ( HXSize_t iPage = 0; iPage < this->GetNPage(); ++ iPage )
            {
                this->GetPage( iPage )->ReSize( 0 );
            }
        }
        return;
    }

    //23 divided by 3 is 7, remainder 2.
    HXSize_t nPage = nLength / maxUnitSize;
    HXLongLong_t remainder = nLength % maxUnitSize;

    HXSize_t additionalPage = 0;
    if ( remainder )
    {
        additionalPage = 1;
    }

    HXSize_t newNPage = nPage + additionalPage;
    this->ResizeNPage( newNPage );

    for ( HXSize_t iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        HXLongLong_t needSize = maxUnitSize;
        if ( iPage == nPage )
        {
            needSize = remainder;
        }
        this->GetPage( iPage )->ReSize( needSize );
    }
}

void DataBook::ResizeNPage( HXSize_t newNPage )
{
    HXSize_t oldNPage = this->GetNPage();

    if ( newNPage <= oldNPage )
    {
        this->Erase( newNPage, oldNPage );
        dataBook->resize( newNPage );
    }
    else
    {
        HXSize_t iPageStart = oldNPage;
        HXSize_t iPageEnd = newNPage;

        for ( HXSize_t iPage = iPageStart; iPage != iPageEnd; ++ iPage )
        {
            dataBook->push_back( new DataPage() );
        }
    }
}

void DataBook::SecureRelativeSpace( HXLongLong_t dataSize )
{
    HXLongLong_t needSize = this->currPos + dataSize;

    this->SecureAbsoluteSpace( needSize );
}

void DataBook::SecureAbsoluteSpace( HXLongLong_t needSize )
{
    //If there is enough space, there is no need to allocate
    //This can cause some std::string problems, and if not ReSize, there may be superfluous characters
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
    for ( HXSize_t iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        this->GetPage( iPage )->MoveToEnd();
    }
}

HXLongLong_t DataBook::GetRemainingSizeOfCurrentPage()
{
    HXLongLong_t remainder = this->currPos % maxUnitSize;
    return maxUnitSize - remainder;
}

void DataBook::ReadFile( std::fstream & file )
{
    //Read the contents of file into DataBook
    //And for DataBook, the process is counter, equivalent to writing

    HXLongLong_t nLength = 0;
    ONEFLOW::HXRead( & file, nLength );

    if ( nLength <= 0 ) return;

    this->SecureAbsoluteSpace( nLength );

    for ( std::streamsize iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        this->GetPage( iPage )->ReadFile( file );
    }
}

void DataBook::WriteFile( std::fstream & file )
{
    HXLongLong_t nLength = this->GetSize();

    //Whether or not nLength is less than zero, you need to write the file
    ONEFLOW::HXWrite( & file, nLength );
    if ( nLength <= 0 )
    {
        return;
    }

    for ( HXSize_t iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        this->GetPage( iPage )->WriteFile( file );
    }
}

void DataBook::ToString( std::string & str )
{
    for ( HXSize_t iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        this->GetPage( iPage )->ToString( str );
    }
}

void DataBook::Append( void * data, HXLongLong_t dataSize )
{
    this->MoveToEnd();
    HXLongLong_t needSize = this->GetSize() + dataSize;
    this->SecureAbsoluteSpace( needSize );

    this->GetCurrentPage()->Write( data, dataSize );
}

void DataBook::Destroy( DataPage * dataPage )
{
    delete dataPage;
}

void DataBook::Erase( HXSize_t startPage, HXSize_t endPage )
{
    for ( HXSize_t iPage = startPage; iPage != endPage; ++ iPage )
    {
        this->Destroy( GetPage( iPage ) );
    }
}

void DataBook::Send( int pid, int tag )
{
    HXLongLong_t nLength = this->GetSize();

    ONEFLOW::HXSend( & nLength, 1, PL_LONG_LONG_INT, pid, tag );

    //It is necessary to judge the zero of data length
    if ( nLength <= 0 ) return;

    for ( HXSize_t iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        this->GetPage( iPage )->Send( pid, tag );
    }
}

void DataBook::Recv( int pid, int tag )
{
    HXLongLong_t nLength = 0;

    ONEFLOW::HXRecv( & nLength, 1, PL_LONG_LONG_INT, pid, tag );

    if ( nLength <= 0 )
    {
        return;
    }

    this->SecureAbsoluteSpace( nLength );

    for ( HXSize_t iPage = 0; iPage < this->GetNPage(); ++ iPage )
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
    HXLongLong_t nLength = this->GetSize();

    HXBcast( & nLength, 1, rootid );

    if ( nLength <= 0 )
    {
        return;
    }

    if ( Parallel::pid != rootid )
    {
        this->SecureAbsoluteSpace( nLength );
    }

    for ( HXSize_t iPage = 0; iPage < this->GetNPage(); ++ iPage )
    {
        this->GetPage( iPage )->Bcast( rootid );
    }
}

void ToDataBook( DataBook * dataBook, std::ostringstream & oss )
{
    if ( ! dataBook ) return;

    dataBook->MoveToBegin();
    dataBook->ReSize( 0 );
    dataBook->Write( & oss );
}

EndNameSpace

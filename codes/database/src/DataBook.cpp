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


#include "DataBook.h"
#include "DataBaseIO.h"
#include "BasicParallel.h"
#include "Parallel.h"
#include <sstream>

BeginNameSpace( ONEFLOW )

DataBook::DataBook(HXOffset_t unitSize)
    : currPos(0), maxUnitSize(unitSize)
{
    if (unitSize <= 0)
        throw std::invalid_argument("DataBook: unitSize must be positive");
}

HXSize_t DataBook::GetPageCount()
{
    return pages.size();
}

DataPage* DataBook::GetCurrentPage()
{
    // Compute on the fly; no need to store currPageId as a member.
    HXSize_t pageId = static_cast<HXSize_t>(this->currPos / maxUnitSize);
    return pages[pageId].get();
}

DataPage * DataBook::GetPage( HXSize_t iPage )
{
    return pages[ iPage ].get();
}

void DataBook::Advance( HXOffset_t offset )
{
    this->currPos += offset;
}

void DataBook::Write( void * data, HXOffset_t dataSize )
{
    if ( dataSize <= 0 ) return;

    this->SecureRelativeSpace( dataSize );

    // Iterative instead of recursive: walk across as many pages as needed
    // in a loop, rather than recursing once per page boundary. This avoids
    // stack depth proportional to (dataSize / maxUnitSize), which becomes
    // a real risk when maxUnitSize is configured small (see
    // DataBook_ManyTinyPages_NoStackOverflow, ~0.54s for just 5000 bytes
    // with unitSize=1 under the old recursive implementation).
    char * cursor = reinterpret_cast<char *>( data );
    HXOffset_t remaining = dataSize;

    while ( remaining > 0 )
    {
        HXOffset_t chunk = this->GetRemainingSizeOfCurrentPage();
        if ( chunk > remaining )
        {
            chunk = remaining;
        }

        this->GetCurrentPage()->Write( cursor, chunk );
        this->Advance( chunk );

        cursor    += chunk;
        remaining -= chunk;
    }
}

void DataBook::Read( void * data, HXOffset_t dataSize )
{
    if ( dataSize <= 0 ) return;

    char * cursor = reinterpret_cast<char *>( data );
    HXOffset_t remaining = dataSize;

    while ( remaining > 0 )
    {
        HXOffset_t chunk = this->GetRemainingSizeOfCurrentPage();
        if ( chunk > remaining )
        {
            chunk = remaining;
        }

        this->GetCurrentPage()->Read( cursor, chunk );
        this->Advance( chunk );

        cursor    += chunk;
        remaining -= chunk;
    }
}

void DataBook::WriteString( std::string & cs )
{
    HXSize_t nLength = cs.length();
    this->Write( & nLength, sizeof( HXSize_t ) );

    // Write the raw bytes directly from the string's own buffer;
    // no manual new[]/delete[], no dependency on a trailing '\0'.
    // cs.data() has been guaranteed contiguous since C++11.
    if ( nLength > 0 )
    {
        this->Write( const_cast<char *>( cs.data() ), nLength );
    }
}

void DataBook::ReadString( std::string & cs )
{
    HXSize_t nLength = 0;
    this->Read( & nLength, sizeof( HXSize_t ) );

    // Resize the string first so it owns the buffer we read into.
    // No manual allocation, and no assumption about a null terminator.
    cs.resize( nLength );
    if ( nLength > 0 )
    {
        this->Read( &cs[0], nLength );
    }
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

HXOffset_t DataBook::GetSize()
{
    HXOffset_t sum = 0;
    for ( int iPage = 0; iPage < this->GetPageCount(); ++ iPage )
    {
        sum += this->GetPage( iPage )->size();
    }
    return sum;
}

// ReSize: align the internal structure with the vector of pages
void DataBook::ReSize(HXOffset_t nLength)
{
    if (nLength < 0) return;          // »ò throw

    if (nLength == 0)
    {
        pages.clear();
        currPos = 0;
        return;
    }

    // Compute how many pages are needed
    // 23 divided by 3 is 7, remainder 2.
    HXSize_t nPage = static_cast<HXSize_t>(nLength / maxUnitSize);
    HXOffset_t remainder = nLength % maxUnitSize;
    HXSize_t newNPage = nPage + (remainder ? 1 : 0);

    SetPageCount(newNPage);       // grows/shrinks the page vector (make_unique when growing)

    // Set the actual size of each page
    for (HXSize_t i = 0; i < newNPage; ++i)
    {
        HXOffset_t pageSize = (i == nPage) ? remainder : maxUnitSize;
        // Edge case: when remainder == 0, the last page is full-sized
        if (i == newNPage - 1 && remainder == 0)
            pageSize = maxUnitSize;
        GetPage(i)->ReSize(static_cast<HXSize_t>(pageSize));
    }

    // Optional: clamp currPos if it exceeds the new size
    if (currPos > nLength)
        currPos = nLength;
}

//void DataBook::SetPageCount(HXSize_t newNPage)
//{
//    HXSize_t oldNPage = pages.size();
//
//    if (newNPage < oldNPage)
//    {
//        // Shrinking: unique_ptr destructor automatically deletes the DataPage.
//        pages.resize(newNPage);
//    }
//    else if (newNPage > oldNPage)
//    {
//        // Growing: must explicitly create real DataPage objects.
//        // vector::resize only default-constructs empty unique_ptrs (nullptr).
//        pages.reserve(newNPage);   // optional, avoids reallocation
//        for (HXSize_t i = oldNPage; i < newNPage; ++i)
//        {
//            pages.push_back(std::make_unique<DataPage>());
//        }
//    }
//    // newNPage == oldNPage ¡ú do nothing
//}

// Sets the number of DataPage objects to exactly newPageCount, growing
// or shrinking as needed. Existing pages below the new count are left
// untouched; existing pages above it are destroyed. Newly created pages
// (when growing) start out empty (size 0) -- callers are responsible for
// sizing their content afterward (see DataBook::ReSize, which calls this
// then sets each page's byte size separately).
void DataBook::SetPageCount( HXSize_t newPageCount )
{
    HXSize_t oldPageCount = pages.size();

    if ( newPageCount < oldPageCount )
    {
        // Shrinking: unique_ptr's destructor automatically deletes the
        // DataPage objects being dropped.
        pages.resize( newPageCount );
    }
    else if ( newPageCount > oldPageCount )
    {
        // Growing: vector::resize() alone only default-constructs empty
        // (nullptr) unique_ptrs, so real DataPage objects must be created
        // explicitly here.
        //
        // Not merely a micro-optimization: without this reserve, growing
        // via repeated emplace_back() could trigger multiple reallocations
        // as capacity doubles, each one moving every existing unique_ptr
        // element to new storage. A single reserve() upfront avoids that
        // entirely, since the target size is already known.
        pages.reserve( newPageCount );
        for ( HXSize_t i = oldPageCount; i < newPageCount; ++ i )
        {
            pages.emplace_back( std::make_unique<DataPage>() );
        }
    }
    // newPageCount == oldPageCount -> nothing to do
}

void DataBook::SecureRelativeSpace( HXOffset_t dataSize )
{
    HXOffset_t needSize = this->currPos + dataSize;

    this->SecureAbsoluteSpace( needSize );
}

void DataBook::SecureAbsoluteSpace( HXOffset_t needSize )
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
    for ( int iPage = 0; iPage < this->GetPageCount(); ++ iPage )
    {
        this->GetPage( iPage )->MoveToBegin();
    }
}

void DataBook::MoveToEnd()
{
    this->currPos = this->GetSize();
    for ( HXSize_t iPage = 0; iPage < this->GetPageCount(); ++ iPage )
    {
        this->GetPage( iPage )->MoveToEnd();
    }
}

HXOffset_t DataBook::GetRemainingSizeOfCurrentPage()
{
    // Offset of the cursor within the current page (0 <= offset < maxUnitSize).
    // e.g. currPos=27, maxUnitSize=10 -> offset=7, meaning we are 7 bytes
    // into page index 2.
    HXOffset_t offsetInPage = this->currPos % maxUnitSize;

    // Bytes left before hitting the end of the current page.
    return maxUnitSize - offsetInPage;
}


void DataBook::ReadFile( std::fstream & file )
{
    //Read the contents of file into DataBook
    //And for DataBook, the process is counter, equivalent to writing

    HXOffset_t nLength = 0;
    ONEFLOW::HXRead( & file, nLength );

    if ( nLength <= 0 ) return;

    this->SecureAbsoluteSpace( nLength );

    for ( std::streamsize iPage = 0; iPage < this->GetPageCount(); ++ iPage )
    {
        this->GetPage( iPage )->ReadFile( file );
    }
}

void DataBook::WriteFile( std::fstream & file )
{
    HXOffset_t nLength = this->GetSize();

    //Whether or not nLength is less than zero, you need to write the file
    ONEFLOW::HXWrite( & file, nLength );
    if ( nLength <= 0 )
    {
        return;
    }

    for ( HXSize_t iPage = 0; iPage < this->GetPageCount(); ++ iPage )
    {
        this->GetPage( iPage )->WriteFile( file );
    }
}

void DataBook::ToString( std::string & str )
{
    for ( HXSize_t iPage = 0; iPage < this->GetPageCount(); ++ iPage )
    {
        this->GetPage( iPage )->ToString( str );
    }
}

void DataBook::Append( void * data, HXOffset_t dataSize )
{
    this->MoveToEnd();
    HXOffset_t needSize = this->GetSize() + dataSize;
    this->SecureAbsoluteSpace( needSize );

    this->GetCurrentPage()->Write( data, dataSize );
}

void DataBook::Send( int pid, int tag )
{
    HXOffset_t nLength = this->GetSize();

    ONEFLOW::HXSend( & nLength, 1, PL_LONG_LONG_INT, pid, tag );

    //It is necessary to judge the zero of data length
    if ( nLength <= 0 ) return;

    for ( HXSize_t iPage = 0; iPage < this->GetPageCount(); ++ iPage )
    {
        this->GetPage( iPage )->Send( pid, tag );
    }
}

void DataBook::Recv( int pid, int tag )
{
    HXOffset_t nLength = 0;

    ONEFLOW::HXRecv( & nLength, 1, PL_LONG_LONG_INT, pid, tag );

    if ( nLength <= 0 )
    {
        return;
    }

    this->SecureAbsoluteSpace( nLength );

    for ( HXSize_t iPage = 0; iPage < this->GetPageCount(); ++ iPage )
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
    HXOffset_t nLength = this->GetSize();

    HXBcast( & nLength, 1, rootid );

    if ( nLength <= 0 )
    {
        return;
    }

    if ( Parallel::pid != rootid )
    {
        this->SecureAbsoluteSpace( nLength );
    }

    for ( HXSize_t iPage = 0; iPage < this->GetPageCount(); ++ iPage )
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

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

#include "BasicParallel.h"
#include <vector>
#include <iostream>
#include <stdexcept>


BeginNameSpace( ONEFLOW )

#ifdef HX_PARALLEL
namespace
{
bool ownsMpiRuntime = false;
}
#endif

void HXFinalize()
{
#ifdef HX_PARALLEL
    if ( !ownsMpiRuntime ) return;

    int finalized = 0;
    const int queryError = MPI_Finalized( &finalized );
    if ( queryError != MPI_SUCCESS )
    {
        throw std::runtime_error( "MPI_Finalized failed during OneFLOW shutdown." );
    }
    if ( !finalized )
    {
        const int finalizeError = MPI_Finalize();
        if ( finalizeError != MPI_SUCCESS )
        {
            throw std::runtime_error( "MPI_Finalize failed during OneFLOW shutdown." );
        }
    }
    ownsMpiRuntime = false;
#endif
}

int HXInit()
{
    int argc = 0;
    char ** argv = nullptr;
    return ONEFLOW::HXInit( argc, &argv );
}

int HXInit( int & argc, char *** argv )
{
#ifdef HX_PARALLEL
    int initialized = 0;
    int error = MPI_Initialized( &initialized );
    if ( error != MPI_SUCCESS ) return error;

    if ( initialized )
    {
        int finalized = 0;
        error = MPI_Finalized( &finalized );
        if ( error != MPI_SUCCESS ) return error;
        return finalized ? MPI_ERR_OTHER : MPI_SUCCESS;
    }

    error = MPI_Init( &argc, argv );
    if ( error == MPI_SUCCESS ) ownsMpiRuntime = true;
    return error;
#else
    (void)argc;
    (void)argv;
    return 0;
#endif
}

int HXRank()
{
    int rank = 0;
#ifdef HX_PARALLEL
    MPI_Comm_rank( MPI_COMM_WORLD, & rank );
#endif
    return rank;
}

int HXSize()
{
    int size = 1;
#ifdef HX_PARALLEL
    MPI_Comm_size( MPI_COMM_WORLD, & size );
#endif
    return size;
}

std::string HXGetProcessorName()
{
    std::string procName = "";
#ifdef HX_PARALLEL
    char cName[ MPI_MAX_PROCESSOR_NAME ];
    int nLength = 0;

    MPI_Get_processor_name( cName, & nLength );
    procName = cName;
#endif
    return procName;
}

void HXSend( const void * data, int size, PL_Datatype dataType, int pid, int tag )
{
#ifdef HX_PARALLEL
    if ( size <= 0 ) return;
    MPI_Send( data, size, dataType, pid, tag, MPI_COMM_WORLD );
#endif
}

void HXRecv( void * data, int size, PL_Datatype dataType, int pid, int tag )
{
#ifdef HX_PARALLEL
    if ( size <= 0 ) return;

    MPI_Status status;
    MPI_Recv( data, size, dataType, pid, tag, MPI_COMM_WORLD, & status );
#endif
}

void HXSendChar( const void * data, int size, int pid, int tag )
{
#ifdef HX_PARALLEL
    //Prevent special circumstances
    if ( size <= 0 ) return;
    MPI_Send( data, size, MPI_CHAR, pid, tag, MPI_COMM_WORLD );
#endif
}

void HXRecvChar( void * data, int size, int pid, int tag )
{
#ifdef HX_PARALLEL
    //Prevent special circumstances
    if ( size <= 0 ) return;

    MPI_Status status;
    MPI_Recv( data, size, MPI_CHAR, pid, tag, MPI_COMM_WORLD, & status );
#endif
}

int HXWait( PL_HXRequest * request )
{
    int errorCode = 0;

#ifdef HX_PARALLEL
    errorCode = MPI_Wait( request, MPI_STATUS_IGNORE );
#endif

    return errorCode;
}

int HXWait( int count, PL_HXRequest * arrayOfRequests )
{
    int errorCode = 0;

    if ( count <= 0 ) return - 1;

#ifdef HX_PARALLEL
    errorCode = MPI_Waitall( count, arrayOfRequests, MPI_STATUSES_IGNORE );
#endif

    return errorCode;
}

void HXSendString( std::string & cs, int pid, int tag )
{
    int nLength = cs.length();

    ONEFLOW::HXSend( & nLength, 1, PL_INT, pid, tag );

    int nLength1 = nLength + 1;
    std::vector<char> data( static_cast<std::size_t>( nLength1 ), '\0' );
    cs.copy( data.data(), nLength );

    ONEFLOW::HXSend( data.data(), nLength1, PL_CHAR, pid, tag );
}

void HXRecvString( std::string & cs, int pid, int tag )
{
    int nLength = 0;
    ONEFLOW::HXRecv( & nLength, 1, PL_INT, pid, tag );

    int nLength1 = nLength + 1;
    std::vector<char> data( static_cast<std::size_t>( nLength1 ), '\0' );

    ONEFLOW::HXRecv( data.data(), nLength1, PL_CHAR, pid, tag );

    cs = data.data();
}

void HXReduceInt( const void * s, void * t, int nElem, PL_Op op )
{
#ifdef HX_PARALLEL
    //MPI_INTEGER error !!!!!!!!!!!!!!!!!!!!!
    //MPI_Allreduce( s, t, nElem, MPI_INTEGER, op, MPI_COMM_WORLD );
    MPI_Allreduce( s, t, nElem, MPI_INT, op, MPI_COMM_WORLD );
#endif

}

void HXReduceReal( const void * s, void * t, int nElem, PL_Op op )
{
#ifdef HX_PARALLEL
    MPI_Allreduce( s, t, nElem, MPI_DOUBLE, op, MPI_COMM_WORLD );
#endif

}


EndNameSpace

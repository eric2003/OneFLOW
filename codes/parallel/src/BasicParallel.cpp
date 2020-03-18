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

#include "BasicParallel.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

void HXFinalize()
{
#ifdef HX_PARALLEL
    int err = MPI_Finalize();
#endif
}

int HXInit()
{
    int err      = 0;
    int argc     = 0;
    char *** argv = 0;
    return ONEFLOW::HXInit( argc, argv );
}

int HXInit( int & argc, char *** argv )
{
    int err = 0;
#ifdef HX_PARALLEL
    err = MPI_Init( & argc, argv );
#endif
    return err;
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
    string procName = "";
#ifdef HX_PARALLEL
    char cName[ MPI_MAX_PROCESSOR_NAME ];
    int nLength = 0;

    MPI_Get_processor_name( cName, & nLength );
    procName = cName;
#endif
    return procName;
}

void HXSend( void * data, int size, PL_Datatype dataType, int pid, int tag )
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

void HXSendChar( void * data, int size, int pid, int tag )
{
#ifdef HX_PARALLEL
    //防止特殊情况
    if ( size <= 0 ) return;
    MPI_Send( data, size, MPI_CHAR, pid, tag, MPI_COMM_WORLD );
#endif
}

void HXRecvChar( void * data, int size, int pid, int tag )
{
#ifdef HX_PARALLEL
    //防止特殊情况
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

void HXSendString( string & cs, int pid, int tag )
{
    int nLength = cs.length();

    ONEFLOW::HXSend( & nLength, 1, PL_INT, pid, tag );

    int nLength1 = nLength + 1;
    char * data = new char[ nLength1 ];

    cs.copy( data, nLength );
    data[ nLength ] = '\0';

    ONEFLOW::HXSend( data, nLength1, PL_CHAR, pid, tag );

    delete[] data;
}

void HXRecvString( string & cs, int pid, int tag )
{
    int nLength = 0;
    ONEFLOW::HXRecv( & nLength, 1, PL_INT, pid, tag );

    int nLength1 = nLength + 1;
    char * data = new char[ nLength1 ];

    ONEFLOW::HXRecv( data, nLength1, PL_CHAR, pid, tag );

    cs = data;

    delete[] data;
}

void HXReduceInt( void * s, void * t, int nElem, PL_Op op )
{
#ifdef HX_PARALLEL
    //MPI_INTEGER error !!!!!!!!!!!!!!!!!!!!!
    //MPI_Allreduce( s, t, nElem, MPI_INTEGER, op, MPI_COMM_WORLD );
    MPI_Allreduce( s, t, nElem, MPI_INT, op, MPI_COMM_WORLD );
#endif

}

void HXReduceReal( void * s, void * t, int nElem, PL_Op op )
{
#ifdef HX_PARALLEL
    MPI_Allreduce( s, t, nElem, MPI_DOUBLE, op, MPI_COMM_WORLD );
#endif

}


EndNameSpace
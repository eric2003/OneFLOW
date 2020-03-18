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

#pragma once

#ifdef HX_PARALLEL
#include "mpi.h"
#endif

#include "Configure.h"
#include <string>
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

#ifdef HX_PARALLEL
    typedef  MPI_Request  PL_HXRequest;
    typedef  MPI_Op       PL_Op;
    typedef  MPI_Datatype PL_Datatype;

    #define PL_REQUEST_NULL    MPI_REQUEST_NULL
    #define PL_MAX             MPI_MAX
    #define PL_MIN             MPI_MIN
    #define PL_SUM             MPI_SUM
    #define PL_CHAR            MPI_CHAR
    #define PL_INT             MPI_INT
    #define PL_LONG_LONG_INT   MPI_LONG_LONG_INT
#else
    typedef  int  PL_HXRequest;
    typedef  int  PL_Op;
    typedef  int  PL_Datatype;

    #define PL_REQUEST_NULL    0
    #define PL_MAX             0
    #define PL_MIN             0
    #define PL_SUM             0
    #define PL_CHAR            0
    #define PL_INT             0
    #define PL_LONG_LONG_INT   0
#endif

int HXInit();
int HXInit( int & argc, char *** argv );
void HXFinalize();

int HXRank();
int HXSize();

std::string HXGetProcessorName();

void HXSend( void * data, int size, PL_Datatype dataType, int pid, int tag = 0 );
void HXRecv( void * data, int size, PL_Datatype dataType, int pid, int tag = 0 );

void HXSendChar( void * data, int size, int pid, int tag = 0 );
void HXRecvChar( void * data, int size, int pid, int tag = 0 );

int HXWait( PL_HXRequest * request );
int HXWait( int count, PL_HXRequest * arrayOfRequests );

void HXSendString( string & cs, int pid, int tag );
void HXRecvString( string & cs, int pid, int tag );

template< typename T >
void HXSmartSend( T * field, int nElement, int pid, int tag )
{
    int bufferSize = nElement * sizeof( T );
    HXSendChar( field, bufferSize, pid, tag );
}

template< typename T >
void HXSmartRecv( T * field, int nElement, int pid, int tag )
{
    int bufferSize = nElement * sizeof( T );
    HXRecvChar( field, bufferSize, pid, tag );
}


void HXReduceInt( void * s, void * t, int nElem, PL_Op op );
void HXReduceReal( void * s, void * t, int nElem, PL_Op op );


EndNameSpace
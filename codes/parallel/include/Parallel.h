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

#include "BasicParallel.h"
#include "DataBaseIO.h"
#include <string>
#include <vector>
using namespace std;

BeginNameSpace( ONEFLOW )

void SetUpParallelEnvironment();

void HXBcast( DataBook * dataBook, int rootid );
void HXBcast( DATA_COMPRESS dataCompression, DATA_DECOMPRESS dataDecompression, int rootid );

template< typename T >
void HXReadBcast( fstream & file, T * field, int nElement, int pid );

template< typename T >
void HXBcast( T * field, int nElement, int pid );

class Parallel
{
public:
    Parallel();
    ~Parallel();
public:
    static void SetDefaultParallelParameter();
    static void SetServerid( int serverid );
    static int GetServerid();
    static void SetPid( int pid );
    static int GetPid();
    static void SetNProc( int nProc );
    static int GetNProc();
    static void TestSayHelloFromEveryProcess();
    static void SetDefaultTag( int defaultTagIn );
    static int GetDefaultTag();
    static void CollectString( string & cs, int rootId, int tag );
    static bool IsServer();
    static int GetFid();
    static void GetSrPid( int zid, int & sPid, int & rPid );
public:
    static int zoneMode;
    static int mode;
    static int pid;
    static int nProc;
    static int serverid;
protected:
    static int defaultTag;
};

template< typename T >
void HXReadBcast( fstream & file, T * field, int nElement, int pid )
{
    if ( nElement <= 0 ) return;

    if ( Parallel::pid == pid )
    {
        ONEFLOW::HXRead( & file, field, nElement );
    }

    if ( Parallel::mode == 0 )
    {
        //server mode
        ONEFLOW::HXBcast( field, nElement, pid );
    }
}

template< typename T >
void HXBcast( T * field, int nElement, int pid )
{
    if ( nElement <= 0 ) return;
#ifdef HX_PARALLEL
    int bufferSize = nElement * sizeof( T );
    MPI_Bcast( field, bufferSize, MPI_CHAR, pid, MPI_COMM_WORLD );
#endif
}

void HXBcastString( string & cs, int pid );

template< typename T >
void HXSwapData( T * field, int nElement, int spid, int rpid, int tag = 0 )
{
    if ( nElement <= 0 ) return;

    if ( spid == rpid ) return;

    if ( Parallel::pid == spid )
    {
        ONEFLOW::HXSmartSend( field, nElement, rpid, tag );
    }
    else if ( Parallel::pid == rpid )
    {
        ONEFLOW::HXSmartRecv( field, nElement, spid, tag );
    }
}

class DataBook;

void HXSwapData( DataBook * dataBook, int spid, int rpid, int tag = 0 );

EndNameSpace
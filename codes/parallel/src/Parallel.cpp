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

#include "Parallel.h"
#include "LogFile.h"
#include "BasicParallel.h"
#include "DataBook.h"
#include "OStream.h"
#include "Zone.h"
#include "ZoneState.h"
#include <iostream>
using namespace std;

BeginNameSpace( ONEFLOW )

ostringstream MyStr;

int Parallel::zoneMode = 0;
int Parallel::mode = 0;
int Parallel::pid;
int Parallel::nProc;
int Parallel::serverid;
int Parallel::defaultTag = 0;

Parallel::Parallel()
{
    ;
}

Parallel::~Parallel()
{
    ;
}

void Parallel::SetDefaultParallelParameter()
{
    ONEFLOW::HXInit();
    Parallel::SetServerid( 0 );
    Parallel::SetPid( ONEFLOW::HXRank() );
    Parallel::SetNProc( ONEFLOW::HXSize() );
}

void Parallel::SetServerid( int serverid )
{
    Parallel::serverid = serverid;
}

int Parallel::GetServerid()
{
    return Parallel::serverid;
}

void Parallel::SetPid( int pid )
{
    Parallel::pid = pid;
}

int Parallel::GetPid()
{
    return Parallel::pid;
}

void Parallel::SetNProc( int nProc )
{
    Parallel::nProc = nProc;
}

int Parallel::GetNProc()
{
    return Parallel::nProc;
}

void Parallel::SetDefaultTag( int defaultTagIn )
{
    Parallel::defaultTag = defaultTagIn;
}

int Parallel::GetDefaultTag()
{
    return Parallel::defaultTag;
}

void Parallel::TestSayHelloFromEveryProcess()
{
    int serverid = Parallel::GetServerid();
    ONEFLOW::StrIO.ClearAll();
    ONEFLOW::StrIO << "Hello from process " << pid << " name = " << ONEFLOW::HXGetProcessorName();
    string cs = ONEFLOW::StrIO.str();
    //cout << cs << "\n";

    Parallel::CollectString( cs, serverid, Parallel::GetDefaultTag() );
}

void Parallel::CollectString( string & cs, int rootId, int tag )
{
    if ( Parallel::GetPid() != rootId )
    {
        ONEFLOW::HXSendString( cs, rootId, tag );
    }
    else
    {
        int nProc = Parallel::GetNProc();
        for ( int pid = 0; pid < nProc; ++ pid )
        {
            if ( pid != rootId )
            {
                ONEFLOW::HXRecvString( cs, pid, tag );
            }
            ONEFLOW::logFile << cs << "\n";
        }
    }
}

bool Parallel::IsServer()
{
    return Parallel::GetPid() == Parallel::GetServerid();
}

int Parallel::GetFid()
{
    int fid;
    if ( Parallel::mode == 0 )
    {
        fid = Parallel::serverid;
    }
    else
    {
        fid = Parallel::pid;
    }
    return fid;
}

void Parallel::GetSrPid( int zid, int & sPid, int & rPid )
{
    if ( Parallel::mode == 0 )
    {
        sPid = Parallel::serverid;
        rPid = ZoneState::pid[ zid ];
    }
    else
    {
        sPid = ZoneState::pid[ zid ];
        rPid = ZoneState::pid[ zid ];
    }
}

void SetUpParallelEnvironment()
{
    Parallel::SetDefaultParallelParameter();
}

void HXBcast( DataBook * dataBook, int rootid )
{
    dataBook->Bcast( rootid );
}

void HXBcast( DATA_COMPRESS dataCompression, DATA_DECOMPRESS dataDecompression, int rootid )
{
    int nProc = Parallel::GetNProc();

    if ( nProc <= 1 ) return;

    DataBook * dataBook = new DataBook();

    if ( Parallel::GetPid() == rootid )
    {
        //Compress data, or store data to dataBook
        dataCompression( dataBook );
    }

    //Pass the dataBook to the required processes
    ONEFLOW::HXBcast( dataBook, rootid );

    if ( Parallel::GetPid() != rootid )
    {
        //Extract the data from dataBook to obtain the required information
        dataDecompression( dataBook );
    }

    delete dataBook;
}

void HXSwapData( DataBook * dataBook, int spid, int rpid, int tag )
{
    if ( spid == rpid ) return;

    if ( Parallel::pid == spid )
    {
        dataBook->Send( rpid, tag );
    }
    else if ( Parallel::pid == rpid )
    {
        dataBook->Recv( spid, tag );
    }
}

void HXBcastString( string & cs, int pid )
{
    int nlen = -1;
    if ( pid == Parallel::pid )
    {
        nlen = cs.length();
    }
    HXBcast( & nlen, 1, pid );

    int nlen1 = nlen + 1;

    char * data = new char[ nlen1 ];

    cs.copy( data, nlen );

    HXBcast( data, nlen, pid );

    if ( pid != Parallel::pid )
    {
        data[ nlen ] = '\0';
        cs = data;
    }

    delete[] data;
}

EndNameSpace
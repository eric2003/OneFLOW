#pragma once
#ifdef HX_PARALLEL
#include "mpi.h"
#endif
#include <vector>
#include "LogFile.h"

class Parallel
{
public:
    static int pid;
    static int nProc;
    static int serverid;
    static int tag;
public:
    static void Init();
    static void Finalize();
public:
    static bool IsServer();
};

void HXSendChar( void * data, int size, int pid, int tag );
void HXRecvChar( void * data, int size, int pid, int tag );

template< typename T >
void HXSendRecvData( T * field, int nElement, int send_pid, int recv_pid, int tag = 0 )
{
    logFile << " send_pid = " << send_pid << "\n";
    logFile << " recv_pid = " << recv_pid << "\n";
    logFile << " Parallel::pid = " << Parallel::pid << "\n";

    if ( send_pid == recv_pid ) return;

    if ( nElement <= 0 ) return;

    int buffer_size = nElement * sizeof( T );

    logFile << "nElement = " << nElement << " buffer_size = " << buffer_size << "\n";

    if ( Parallel::pid == send_pid )
    {
        HXSendChar( field, buffer_size, recv_pid, tag );
    }
    else if ( Parallel::pid == recv_pid )
    {
        HXRecvChar( field, buffer_size, send_pid, tag );
    }
}
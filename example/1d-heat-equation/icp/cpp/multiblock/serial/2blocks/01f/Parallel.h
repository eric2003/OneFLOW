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

void HXSendChar( void * data, int size, int pid, int tag = 0 );
void HXRecvChar( void * data, int size, int pid, int tag = 0 );
void HXSendRecvString( std::string & str, int send_pid, int recv_pid, int tag = 0 );
void HXSendString( std::string const & str, int recv_pid, int tag = 0 );
void HXRecvString( std::string & str, int send_pid, int tag = 0 );

template< typename T >
void HXSendData( T * field, int nElement, int recv_pid, int tag  );

template< typename T >
void HXRecvData( T * field, int nElement, int send_pid, int tag );

template< typename T >
void HXSendData( T * field, int nElement, int recv_pid, int tag = 0 )
{
    if ( nElement <= 0 ) return;

    int buffer_size = nElement * sizeof( T );

    HXSendChar( field, buffer_size, recv_pid, tag );
}

template< typename T >
void HXRecvData( T * field, int nElement, int send_pid, int tag = 0 )
{
    if ( nElement <= 0 ) return;

    int buffer_size = nElement * sizeof( T );

    HXRecvChar( field, buffer_size, send_pid, tag );
}

template< typename T >
void HXSendRecvData( T * field, int nElement, int send_pid, int recv_pid, int tag = 0 )
{
    if ( send_pid == recv_pid ) return;

    if ( nElement <= 0 ) return;

    int buffer_size = nElement * sizeof( T );

    if ( Parallel::pid == send_pid )
    {
        HXSendChar( field, buffer_size, recv_pid, tag );
    }
    else if ( Parallel::pid == recv_pid )
    {
        HXRecvChar( field, buffer_size, send_pid, tag );
    }
}

template< typename T >
void HXBcastData( T * field, int nElement, int send_pid )
{
    if ( nElement <= 0 ) return;
    int buffer_size = nElement * sizeof( T );
#ifdef HX_PARALLEL
    MPI_Bcast( field, buffer_size, MPI_CHAR, send_pid, MPI_COMM_WORLD );
#endif
}

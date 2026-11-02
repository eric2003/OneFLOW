#include "Parallel.h"
#include <iostream>

int Parallel::pid = 0;
int Parallel::nProc = 1;
int Parallel::serverid = 0;
int Parallel::tag = 0;

void Parallel::Init()
{
#ifdef HX_PARALLEL
    int argc     = 0;
    char ** argv = 0;
    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &Parallel::pid );
    MPI_Comm_size( MPI_COMM_WORLD, &Parallel::nProc );
    int len = -1;
    char version[ MPI_MAX_LIBRARY_VERSION_STRING ];
    MPI_Get_library_version( version, &len );
    std::cout << "Hello, world!  I am " << Parallel::pid << " of " << Parallel::nProc
        << "(" << version << ", " << len << ")" << std::endl;
#endif
}

void Parallel::Finalize()
{
#ifdef HX_PARALLEL
    MPI_Finalize();
#endif
}

bool Parallel::IsServer()
{
    return Parallel::pid == Parallel::serverid;
}

void HXSendChar( void * data, int size, int pid, int tag )
{
#ifdef HX_PARALLEL
    if ( size <= 0 ) return;
    MPI_Send( data, size, MPI_CHAR, pid, tag, MPI_COMM_WORLD );
#endif
}

void HXRecvChar( void * data, int size, int pid, int tag )
{
#ifdef HX_PARALLEL
    if ( size <= 0 ) return;

    MPI_Status status;
    MPI_Recv( data, size, MPI_CHAR, pid, tag, MPI_COMM_WORLD, & status );
#endif
}

void HXSendString( std::string const & str, int recv_pid, int tag )
{
#ifdef HX_PARALLEL
    unsigned len = str.size();
    MPI_Send( &len, 1, MPI_UNSIGNED, recv_pid, tag, MPI_COMM_WORLD );
    if ( len == 0 ) return;
    MPI_Send( str.data(), len, MPI_CHAR, recv_pid, tag, MPI_COMM_WORLD );
#endif
}

void HXRecvString( std::string & str, int send_pid, int tag )
{
#ifdef HX_PARALLEL
    unsigned len;
    MPI_Status status;
    MPI_Recv( &len, 1, MPI_UNSIGNED, send_pid, tag, MPI_COMM_WORLD, &status );
    if ( len == 0 ) return;
    str.resize( len );
    MPI_Recv( str.data(), len, MPI_CHAR, send_pid, tag, MPI_COMM_WORLD, &status );
#endif
}

void HXSendRecvString( std::string & str, int send_pid, int recv_pid, int tag )
{
    if ( send_pid == recv_pid ) return;
    if ( Parallel::pid == send_pid )
    {
        HXSendString( str, recv_pid, tag );
    }
    else if ( Parallel::pid == recv_pid )
    {
        HXRecvString( str, send_pid, tag );
    }
}


//void HXSendRecvString( std::string & str, int send_pid, int recv_pid, int tag )
//{
//    if ( send_pid == recv_pid ) return;
//    int nsize = -1;
//    if ( Parallel::pid == send_pid )
//    {
//        nsize = str.size();
//    }
//    HXSendRecvData( &nsize, 1, send_pid, recv_pid, tag );
//    if ( Parallel::pid == recv_pid )
//    {
//        str.resize( nsize );
//    }
//    HXSendRecvData( str.data(), nsize, send_pid, recv_pid, tag );
//}
//

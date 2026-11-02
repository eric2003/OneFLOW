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


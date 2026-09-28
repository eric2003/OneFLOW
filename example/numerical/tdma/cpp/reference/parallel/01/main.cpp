#include "mpi.h"
import std;

int main( int argc, char ** argv )
{
    //MPI initialization and both memory allocation and generation of the matrix A
    int i, j, k, size, index;
    int index1, index2;
    int mynode, totalnodes;
    double alpha, gamma;
    const int numrows = 5;

    MPI_Status status;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &totalnodes );
    MPI_Comm_rank( MPI_COMM_WORLD, &mynode );
    size = (int)std::pow( 2, std::log2( totalnodes + 1 ) + 1 ) - 1;
    std::print( "size={}", size );
    double ** A = new double * [ numrows ];
    for ( i = 0; i < numrows; i++ )
    {
        A[ i ] = new double[ size + 1 ];
        for ( j = 0; j < size + 1; j++ )
        {
            A[ i ][ j ] = 0.0;
        }
    }
    if ( mynode == 0 )
    {
        A[ 0 ][ 0 ] = -2.0; A[ 0 ][ 1 ] = 1.0;
        A[ 1 ][ 0 ] = 1.0; A[ 1 ][ 1 ] = -2.0; A[ 1 ][ 2 ] = 1.0;
        A[ 2 ][ 1 ] = 1.0; A[ 2 ][ 2 ] = -2.0; A[ 2 ][ 3 ] = 1.0;
    }
    else if ( mynode == ( totalnodes - 1 ) )
    {
        index = 2 * mynode;
        A[ 0 ][ index - 1 ] = 1.0; A[ 0 ][ index ] = -2.0;
        A[ 0 ][ index + 1 ] = 1.0;
        index = 2 * mynode + 1;
        A[ 1 ][ index - 1 ] = 1.0; A[ 1 ][ index ] = -2.0;
        A[ 1 ][ index + 1 ] = 1.0;
        A[ 2 ][ size - 2 ] = 1.0; A[ 2 ][ size - 1 ] = -2.0;
    }
    else
    {
        for ( i = 0; i < 3; i++ )
        {
            index = i + 2 * mynode;
            A[ i ][ index - 1 ] = 1.0;
            A[ i ][ index ] = -2.0;
            A[ i ][ index + 1 ] = 1.0;
        }
    }
    for ( i = 0; i < 3; i++ )
    {
        A[ i ][ size ] = 2 * mynode + i;
    }
    int numactivep = totalnodes;
    int * activep = new int[ totalnodes ];
    for ( j = 0; j < numactivep; j++ )
    {
        activep[ j ] = j;
    }

    for ( j = 0; j < size + 1; j++ )
    {
        A[ 3 ][ j ] = A[ 0 ][ j ];
        A[ 4 ][ j ] = A[ 2 ][ j ];
    }

    /* Part 3 */
    double * x = new double[ totalnodes ];

    delete[] activep;
    for ( i = 0; i < numrows; i++ )
    {
        delete[] A[ i ];
    }
    delete[] A;
    delete[] x;
    MPI_Finalize();
    return 0;
}

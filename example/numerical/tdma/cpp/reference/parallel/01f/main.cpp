#include "mpi.h"
#include "LogFile.h"
import std;

void PrintToScreen( double ** A, int numrows, int N, const std::string & name = "matrix" )
{
    std::print( "{} = \n", name );
    for ( int i = 0; i < numrows; ++ i )
    {
        for ( int j = 0; j < N; ++ j )
        {
            std::print( "{} ", A[ i ][ j ] );
        }
        std::println();
    }
}

void Print( double ** A, int numrows, int N, const std::string & name = "matrix" )
{
    std::print( "{} = \n", name );
    logFile << name << " = \n";
    for ( int i = 0; i < numrows; ++ i )
    {
        for ( int j = 0; j < N; ++ j )
        {
            logFile << A[ i ][ j ] << " ";
        }
        logFile << "\n";
    }
}

int main( int argc, char ** argv )
{
    //MPI initialization and both memory allocation and generation of the matrix A
    int index;
    int index1, index2;
    double alpha, gamma;
    const int numrows = 5;

    MPI_Status status;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &Parallel::nprocs );
    MPI_Comm_rank( MPI_COMM_WORLD, &Parallel::pid );

    int N = (int)std::pow( 2, std::log2( Parallel::nprocs + 1 ) + 1 ) - 1;
    std::print( "N={}", N );
    double ** A = new double * [ numrows ];
    for ( int i = 0; i < numrows; ++ i )
    {
        A[ i ] = new double[ N + 1 ];
        for ( int j = 0; j < N + 1; ++ j )
        {
            A[ i ][ j ] = 0.0;
        }
    }
    std::print( "Parallel::pid={}, Parallel::nprocs={}, numrows={}, N={}\n", Parallel::pid, Parallel::nprocs, numrows, N );
    Print( A, numrows, N + 1, "Matrix A 0" );
    if ( Parallel::pid == 0 )
    {
        A[ 0 ][ 0 ] = -2.0;
        A[ 0 ][ 1 ] = 1.0;
        A[ 1 ][ 0 ] = 1.0;
        A[ 1 ][ 1 ] = -2.0;
        A[ 1 ][ 2 ] = 1.0;
        A[ 2 ][ 1 ] = 1.0; 
        A[ 2 ][ 2 ] = -2.0; 
        A[ 2 ][ 3 ] = 1.0;
    }
    else if ( Parallel::pid == ( Parallel::nprocs - 1 ) )
    {
        index = 2 * Parallel::pid;
        A[ 0 ][ index - 1 ] = 1.0; 
        A[ 0 ][ index     ] = -2.0;
        A[ 0 ][ index + 1 ] = 1.0;

        index = 2 * Parallel::pid + 1;
        A[ 1 ][ index - 1 ] = 1.0; 
        A[ 1 ][ index     ] = -2.0;
        A[ 1 ][ index + 1 ] = 1.0;
        A[ 2 ][ N - 2 ] = 1.0;
        A[ 2 ][ N - 1 ] = -2.0;
    }
    else
    {
        for ( int i = 0; i < 3; ++ i )
        {
            index = 2 * Parallel::pid + i;
            A[ i ][ index - 1 ] = 1.0;
            A[ i ][ index     ] = -2.0;
            A[ i ][ index + 1 ] = 1.0;
        }
    }
    for ( int i = 0; i < 3; ++ i )
    {
        A[ i ][ N ] = 2 * Parallel::pid + i;
    }
    int numactivep = Parallel::nprocs;
    int * activep = new int[ Parallel::nprocs ];
    for ( int j = 0; j < numactivep; ++ j )
    {
        activep[ j ] = j;
    }

    for ( int j = 0; j < N + 1; ++ j )
    {
        A[ 3 ][ j ] = A[ 0 ][ j ];
        A[ 4 ][ j ] = A[ 2 ][ j ];
    }
    Print( A, numrows, N + 1, "Matrix A 1" );

    /* Part 2 */
    //Cyclic reduction

    for ( int i = 0; i < log2( N + 1 ) - 1; ++ i )
    {
        for ( int j = 0; j < numactivep; ++ j )
        {
            if ( Parallel::pid == activep[ j ] )
            {
                index1 = 2 * Parallel::pid + 1 - pow( 2, i );
                index2 = 2 * Parallel::pid + 1 + pow( 2, i );
                alpha = A[ 1 ][ index1 ] / A[ 3 ][ index1 ];
                gamma = A[ 1 ][ index2 ] / A[ 4 ][ index2 ];
                for ( int k = 0; k < N + 1; ++ k )
                {
                    A[ 1 ][ k ] -= ( alpha * A[ 3 ][ k ] + gamma * A[ 4 ][ k ] );
                }
                if ( numactivep > 1 )
                {
                    if ( j == 0 )
                    {
                        MPI_Send( A[ 1 ], N + 1, MPI_DOUBLE, activep[ 1 ], 0, MPI_COMM_WORLD );
                    }
                    else if ( j == numactivep - 1 )
                    {
                        MPI_Send( A[ 1 ], N + 1, MPI_DOUBLE, activep[ numactivep - 2 ], 1, MPI_COMM_WORLD );
                    }
                    else if ( j % 2 == 0 )
                    {
                        MPI_Send( A[ 1 ], N + 1, MPI_DOUBLE, activep[ j - 1 ], 1, MPI_COMM_WORLD );
                        MPI_Send( A[ 1 ], N + 1, MPI_DOUBLE, activep[ j + 1 ], 0, MPI_COMM_WORLD );
                    }
                    else
                    {
                        MPI_Recv( A[ 3 ], N + 1, MPI_DOUBLE, activep[ j - 1 ], 0, MPI_COMM_WORLD, &status );
                        MPI_Recv( A[ 4 ], N + 1, MPI_DOUBLE, activep[ j + 1 ], 1, MPI_COMM_WORLD, &status );
                    }
                }
            }
        }
        numactivep = 0;
        for ( int j = activep[ 1 ]; j < Parallel::nprocs; j = j + pow( 2, i + 1 ) )
        {
            activep[ numactivep ++ ] = j;
        }
    }

    /* Part 3 */
    //Part 3 - Back substitution
    double * x = new double[ Parallel::nprocs ];
    for ( int j = 0; j < Parallel::nprocs; ++ j )
    {
        x[ j ] = 0.0;
    }
    if ( Parallel::pid == activep[ 0 ] )
    {
        x[ Parallel::pid ] = A[ 1 ][ N ] / A[ 1 ][ ( N - 1 ) / 2 ];
    }
    double tmp;
    for ( int i = log2( N + 1 ) - 3; i >= 0; -- i )
    {
        tmp = x[ Parallel::pid ];
        MPI_Allgather( &tmp, 1, MPI_DOUBLE, x, 1, MPI_DOUBLE, MPI_COMM_WORLD );
        numactivep = 0;
        for ( int j = activep[ 0 ] - pow( 2, i ); j < Parallel::nprocs; j = j + pow( 2, i + 1 ) )
        {
            activep[ numactivep ++ ] = j;
        }
        for ( int j = 0; j < numactivep; ++ j )
        {
            if ( Parallel::pid == activep[ j ] )
            {
                x[ Parallel::pid ] = A[ 1 ][ N ];
                for ( int k = 0; k < Parallel::nprocs; ++ k )
                {
                    if ( k != Parallel::pid )
                    {
                        x[ Parallel::pid ] -= A[ 1 ][ 2 * k + 1 ] * x[ k ];
                    }
                }
                x[ Parallel::pid ] = x[ Parallel::pid ] / A[ 1 ][ 2 * Parallel::pid + 1 ];
            }
        }
    }

    tmp = x[ Parallel::pid ];
    MPI_Allgather( &tmp, 1, MPI_DOUBLE, x, 1, MPI_DOUBLE, MPI_COMM_WORLD );

    /* Part 4 */
    //Part 4 - Solving for odd rows

    for ( int k = 0; k < Parallel::nprocs; ++ k )
    {
        A[ 0 ][ N ] -= A[ 0 ][ 2 * k + 1 ] * x[ k ];
        A[ 2 ][ N ] -= A[ 2 ][ 2 * k + 1 ] * x[ k ];
    }
    A[ 0 ][ N ] = A[ 0 ][ N ] / A[ 0 ][ 2 * Parallel::pid ];
    A[ 1 ][ N ] = x[ Parallel::pid ];
    A[ 2 ][ N ] = A[ 2 ][ N ] / A[ 2 ][ 2 * Parallel::pid + 2 ];

    std::print( "Parallel::nprocs={}\n", Parallel::nprocs );

    if ( Parallel::pid == 0 )
    {
        std::print( "Parallel::pid = {}, x = ", Parallel::pid );
        for ( int j = 0; j < Parallel::nprocs; ++ j )
        {
            std::print( "{} ", x[ j ] );
        }
        std::println();
    }

    Print( A, numrows, N + 1, "Matrix A Final" );


    delete[] activep;
    for ( int i = 0; i < numrows; ++ i )
    {
        delete[] A[ i ];
    }
    delete[] A;
    delete[] x;
    MPI_Finalize();

    return 0;
}

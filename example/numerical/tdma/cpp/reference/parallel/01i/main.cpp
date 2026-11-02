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

void Print( std::vector<std::vector<double>> & a, const std::string & name = "matrix" )
{
    int NI = a.size();
    logFile << name << " = \n";
    for ( int i = 0; i < NI; ++ i )
    {
        int NJ = a[ i ].size();
        for ( int j = 0; j < NJ; ++ j )
        {
            logFile << std::format( "{:5.2f} ", a[ i ][ j ] );
        }
        logFile << "\n";
    }
    logFile << "\n";
}

template<typename T>
void Print( std::vector<T> & x, const std::string & name = "vector" )
{
    logFile << name << " = \n";
    for ( auto v: x )
    {
        logFile << std::format( "{:4} ", v );
    }
    logFile << "\n";
}


int main( int argc, char ** argv )
{
    //Part 1 - MPI initialization
    const int numrows = 5;

    MPI_Status status;
    MPI_Init( &argc, &argv );
    MPI_Comm_size( MPI_COMM_WORLD, &Parallel::nprocs );
    MPI_Comm_rank( MPI_COMM_WORLD, &Parallel::pid );

    int N = (int)std::pow( 2, std::log2( Parallel::nprocs + 1 ) + 1 ) - 1;
    std::print( "N={}\n", N );
    std::print( "std::log2( Parallel::nprocs + 1 ) + 1 = {}\n", std::log2( Parallel::nprocs + 1 ) + 1 );
    std::vector<std::vector<double>> A( numrows );
    for ( int i = 0; i < A.size(); ++ i )
    {
        A[ i ].resize( N + 1, 0.0 );
    }
    std::print( "Parallel::pid={}, Parallel::nprocs={}, numrows={}, N={}\n", Parallel::pid, Parallel::nprocs, numrows, N );
    Print( A, "Matrix A 0" );
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
        int index = 2 * Parallel::pid;
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
            int index = 2 * Parallel::pid + i;
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
    std::vector<int> activep( Parallel::nprocs );
    for ( int j = 0; j < numactivep; ++ j )
    {
        activep[ j ] = j;
    }

    for ( int j = 0; j < N + 1; ++ j )
    {
        A[ 3 ][ j ] = A[ 0 ][ j ];
        A[ 4 ][ j ] = A[ 2 ][ j ];
    }
    Print( A, "Matrix A 1" );

    //Part 2 - Cyclic reduction

    for ( int level = 0; level < std::log2( N + 1 ) - 1; ++ level )
    {
        logFile << "Cyclic reduction level = " << level << "\n";
        int offset = std::pow( 2, level );
        int step = std::pow( 2, level + 1 );
        logFile << "offset = " << offset << " step = " << step << "\n";
        logFile << "numactivep = " << numactivep  << "\n";
        for ( int j = 0; j < numactivep; ++ j )
        {
            logFile << "j = " << j << " numactivep = " << numactivep  << "\n";
            if ( Parallel::pid == activep[ j ] )
            {
                logFile << " Parallel::pid = " << Parallel::pid << " activep[ j ] = " << activep[ j ] << "\n";

                int index1 = 2 * Parallel::pid + 1 - offset;
                int index2 = 2 * Parallel::pid + 1 + offset;
                logFile << " index1 = " << index1 << " index2 = " << index2 << "\n";

                double alpha = A[ 1 ][ index1 ] / A[ 3 ][ index1 ];
                double gamma = A[ 1 ][ index2 ] / A[ 4 ][ index2 ];
                for ( int k = 0; k < N + 1; ++ k )
                {
                    A[ 1 ][ k ] -= ( alpha * A[ 3 ][ k ] + gamma * A[ 4 ][ k ] );
                }
                if ( numactivep > 1 )
                {
                    if ( j == 0 )
                    {
                        MPI_Send( A[ 1 ].data(), N + 1, MPI_DOUBLE, activep[ 1 ], 0, MPI_COMM_WORLD );
                    }
                    else if ( j == numactivep - 1 )
                    {
                        MPI_Send( A[ 1 ].data(), N + 1, MPI_DOUBLE, activep[ numactivep - 2 ], 1, MPI_COMM_WORLD );
                    }
                    else if ( j % 2 == 0 )
                    {
                        MPI_Send( A[ 1 ].data(), N + 1, MPI_DOUBLE, activep[ j - 1 ], 1, MPI_COMM_WORLD );
                        MPI_Send( A[ 1 ].data(), N + 1, MPI_DOUBLE, activep[ j + 1 ], 0, MPI_COMM_WORLD );
                    }
                    else
                    {
                        MPI_Recv( A[ 3 ].data(), N + 1, MPI_DOUBLE, activep[ j - 1 ], 0, MPI_COMM_WORLD, &status );
                        MPI_Recv( A[ 4 ].data(), N + 1, MPI_DOUBLE, activep[ j + 1 ], 1, MPI_COMM_WORLD, &status );
                    }
                }
            }
        }

        logFile << "activep[ 1 ] = " << activep[ 1 ] << " numactivep = " << numactivep;
        logFile << " Parallel::nprocs = " << Parallel::nprocs << "\n";
        logFile << "Cyclic reduction level = " << level << "\n";
        Print( activep, "Cyclic reduction vector activep 000" );

        numactivep = 0;
        for ( int j = activep[ 1 ]; j < Parallel::nprocs; j = j + step )
        {
            logFile << "j = " << j << " numactivep = " << numactivep << "\n";
            activep[ numactivep ++ ] = j;
        }
        logFile << "Final numactivep = " << numactivep << "\n";
        Print( activep, "Cyclic reduction vector activep 111" );
    }

    Print( A, "Matrix A 2" );

    //Part 3 - Back substitution
    logFile << "Back substitution\n";
    std::vector<double> x( Parallel::nprocs, 0.0 );
    logFile << "Parallel::pid = " << Parallel::pid << " activep[ 0 ] = " << activep[ 0 ] << "\n";
    if ( Parallel::pid == activep[ 0 ] )
    {
        logFile << "N = " << N << ";( N - 1 ) / 2 = " << ( N - 1 ) / 2 << "\n";
        x[ Parallel::pid ] = A[ 1 ][ N ] / A[ 1 ][ ( N - 1 ) / 2 ];
    }
    double tmp;
    std::print( "std::log2( N + 1 ) = {}\n", std::log2( N + 1 ) );
    std::print( "std::log2( N + 1 ) - 3 = {}\n", std::log2( N + 1 ) - 3 );
    std::print( "std::log2( N + 1 ) - 2 = {}\n", std::log2( N + 1 ) - 2 );
    Print( x, "vector x 000" );
    for ( int level = std::log2( N + 1 ) - 3; level >= 0; -- level )
    {
        logFile << "Back substitution level = " << level << "\n";
        int offset = std::pow( 2, level );
        int step = std::pow( 2, level + 1 );

        logFile << "offset = " << offset << " step = " << step << "\n";

        tmp = x[ Parallel::pid ];
        MPI_Allgather( &tmp, 1, MPI_DOUBLE, x.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD );
        Print( x, "vector x 111" );
        logFile << "activep[ 0 ] = " << activep[ 0 ] << " activep[ 0 ] - offset = " << activep[ 0 ] - offset;
        logFile << " Parallel::nprocs = " << Parallel::nprocs << "\n";
        Print( activep, "vector activep 000" );
        numactivep = 0;
        for ( int j = activep[ 0 ] - offset; j < Parallel::nprocs; j = j + step )
        {
            logFile << "j = " << j << " numactivep = " << numactivep << "\n";
            logFile << "000 activep[" << numactivep << "] = " << activep[ numactivep ] << "\n";
            activep[ numactivep ++ ] = j;
            logFile << "111 activep[" << numactivep << "] = " << activep[ numactivep ] << "\n";
        }
        logFile << "numactivep = " << numactivep << "\n";
        Print( activep, "vector activep 111" );

        for ( int j = 0; j < numactivep; ++ j )
        {
            logFile << "j = " << j << " numactivep = " << numactivep;
            logFile << " activep[" << j << "] = " << activep[ j ] << " Parallel::pid = " << Parallel::pid << "\n";
                
            if ( Parallel::pid == activep[ j ] )
            {
                logFile << "000 x[" << Parallel::pid << "] = " << x[ Parallel::pid ] << "\n";
                x[ Parallel::pid ] = A[ 1 ][ N ];
                logFile << "111 x[" << Parallel::pid << "] = " << x[ Parallel::pid ] << "\n";
                for ( int k = 0; k < Parallel::nprocs; ++ k )
                {
                    if ( k != Parallel::pid )
                    {
                        x[ Parallel::pid ] -= A[ 1 ][ 2 * k + 1 ] * x[ k ];
                    }
                }
                logFile << "222 x[" << Parallel::pid << "] = " << x[ Parallel::pid ] << "\n";
                logFile << "Parallel::pid = "<< Parallel::pid << "; 2 * Parallel::pid + 1 = " << 2 * Parallel::pid + 1 << "\n";
                logFile << "A[ 1 ][" << 2 * Parallel::pid + 1 << "] = "<< A[ 1 ][ 2 * Parallel::pid + 1 ] << "\n";
                
                x[ Parallel::pid ] = x[ Parallel::pid ] / A[ 1 ][ 2 * Parallel::pid + 1 ];
                logFile << "333 x[" << Parallel::pid << "] = " << x[ Parallel::pid ] << "\n";
            }
        }
        Print( x, "vector x 222" );
    }

    tmp = x[ Parallel::pid ];
    MPI_Allgather( &tmp, 1, MPI_DOUBLE, x.data(), 1, MPI_DOUBLE, MPI_COMM_WORLD );

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
    Print( x, "vector x" );

    if ( Parallel::pid == 0 )
    {
        std::print( "Parallel::pid = {}, x = ", Parallel::pid );
        for ( int j = 0; j < Parallel::nprocs; ++ j )
        {
            std::print( "{} ", x[ j ] );
        }
        std::println();
    }

    Print( A, "Matrix A Final" );

    MPI_Finalize();

    return 0;
}

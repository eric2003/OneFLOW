import std;

void Print( std::vector<double> & x, const std::string & name = "vector" )
{
    std::print( "{} = ", name );
    for ( auto v: x )
    {
        std::print( "{} ", v );
    }
    std::println();
}

void Print( std::vector<std::vector<double>> & a, const std::string & name = "matrix" )
{
    int N = a.size();
    for ( int i = 0; i < N; ++ i )
    {
        for ( int j = 0; j < N; ++ j )
        {
            std::print( "{:10f} ", a[ i ][ j ] );
        }
        std::println();
    }
    std::println();
}

void MatrixMultiply( std::vector<std::vector<double>> & a, std::vector<double> & x, std::vector<double> & y )
{
    int N = a.size();
    for ( int i = 0; i < N; ++ i )
    {
        y[ i ] = 0.0;
        for ( int j = 0; j < N; ++ j )
        {
            y[ i ] += a[ i ][ j ] * x[ j ];
        }
    }
}


int main( int argc, char ** argv )
{
    //Memory allocation and generation of matrix
    const int size = 3;
    //const int size = 15;
    std::vector<double> x( size );
    for ( int i = 0; i < size; ++ i )
    {
        x[ i ] = 0.0;
    }
    std::vector<double> F( size );
    std::vector<std::vector<double>> A( size );
    for ( int i = 0; i < size; ++ i )
    {
        A[ i ].resize( size );
        for ( int j = 0; j < size; ++ j )
        {
            A[ i ][ j ] = 0.0;
        }
        F[ i ] = i;
    }
    A[ 0 ][ 0 ] = -2.0;
    A[ 0 ][ 1 ] = 1.0;
    A[ size - 1 ][ size - 2 ] = 1.0;
    A[ size - 1 ][ size - 1 ] = -2.0;
    for ( int i = 1; i < size - 1; ++ i )
    {
        A[ i ][ i ] = -2.0;
        A[ i ][ i - 1 ] = 1.0;
        A[ i ][ i + 1 ] = 1.0;
    }

    Print( A, "matrix A" );
    Print( F, "vector F" );

    std::vector<std::vector<double>> AA = A;
    std::vector<double> FF = F;

    //Cyclic reduction
    for ( int level = 0; level < std::log2( size + 1 ) - 1; ++ level )
    {
        std::print( "Cyclic reduction level = {}\n", level );
        int joffset = std::pow( 2, level );
        int jstep = std::pow( 2, level + 1 );
        int jstart = jstep - 1;
        for ( int j = jstart; j < size; j = j + jstep )
        {
            int in = j - joffset;
            int ip = j + joffset;
            double alpha = A[ j ][ in ] / A[ in ][ in ];
            double gamma = A[ j ][ ip ] / A[ ip ][ ip ];
            std::print( "in,j,ip = ({},{},{})\n", in, j, ip );
            std::print( "A[{}][{}]={}\n", j, in, A[ j ][ in ] );
            std::print( "A[{}][{}]={}\n", j, ip, A[ j ][ ip ] );
            std::print( "A[{}][{}]={}\n", in, in, A[ in ][ in ] );
            std::print( "A[{}][{}]={}\n", ip, ip, A[ ip ][ ip ] );
            std::print( "alpha = {}, gamma = {}\n", alpha, gamma );
            for ( int k = 0; k < size; ++ k )
            {
                A[ j ][ k ] -= ( alpha * A[ in ][ k ] + gamma * A[ ip ][ k ] );
                std::print( "j,k=({},{}) ", j, k );
                std::print( "in,k=({},{}) ", in, k );
                std::print( "ip,k=({},{}) \n", ip, k );
            }
            std::println();
            F[ j ] -= ( alpha * F[ in ] + gamma * F[ ip ] );
        }
    }

    Print( A, "matrix A after cyclic reduction" );
    Print( F, "vector F after cyclic reduction" );
    //Back substitution
    int index = ( size - 1 ) / 2;
    x[ index ] = F[ index ] / A[ index ][ index ];
    std::print( "F[{}] = {}, A[{}][{}]={}\n", index, F[ index ], index, index, A[ index ][ index ] );
    std::print( "x[{}] = {}\n", index, x[ index ] );
    for ( int level = std::log2( size + 1 ) - 2; level >= 0; -- level )
    {
        std::print( "Back substitution level = {}\n", level );
        int joffset = std::pow( 2, level );
        int jstep = std::pow( 2, level + 1 );
        int jstart = jstep - 1;
        for ( int j = jstart; j < size; j = j + jstep )
        {
            int in = j - joffset;
            int ip = j + joffset;
            std::print( "in,j,ip = ({},{},{})\n", in, j, ip );
            x[ in ] = F[ in ];
            x[ ip ] = F[ ip ];
            for ( int k = 0; k < size; ++ k )
            {
                if ( k != in )
                {
                    x[ in ] -= A[ in ][ k ] * x[ k ];
                }
                if ( k != ip )
                {
                    x[ ip ] -= A[ ip ][ k ] * x[ k ];
                }
            }
            x[ in ] = x[ in ] / A[ in ][ in ];
            x[ ip ] = x[ ip ] / A[ ip ][ ip ];
            std::print( "x[{}] = {}, x[{}] = {}\n", in, x[ in ], ip, x[ ip ] );
        }
    }

    for ( int i = 0; i < size; ++ i )
    {
        std::println( "{}", x[i] );
    }

    std::vector<double> y( F.size() );
    MatrixMultiply( AA, x, y );
    Print( y, "y" );
    Print( FF, "FF" );

    return 0;
}

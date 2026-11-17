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
    //const int size = 15;
    const int size = 3;
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
    for ( int i = 0; i < std::log2( size + 1 ) - 1; ++ i )
    {
        std::print( "Cyclic reduction i = {}\n", i );
        for ( int j = std::pow( 2, i + 1 ) - 1; j < size; j = j + std::pow( 2, i + 1 ) )
        {
            int offset = std::pow( 2, i );
            int index1 = j - offset;
            int index2 = j + offset;
            double alpha = A[ j ][ index1 ] / A[ index1 ][ index1 ];
            double gamma = A[ j ][ index2 ] / A[ index2 ][ index2 ];
            std::print( "index1,j,index2 = ({},{},{})\n", index1, j, index2 );
            std::print( "A[{}][{}]={}\n", j, index1, A[ j ][ index1 ] );
            std::print( "A[{}][{}]={}\n", j, index2, A[ j ][ index2 ] );
            std::print( "A[{}][{}]={}\n", index1, index1, A[ index1 ][ index1 ] );
            std::print( "A[{}][{}]={}\n", index2, index2, A[ index2 ][ index2 ] );
            std::print( "alpha = {}, gamma = {}\n", alpha, gamma );
            for ( int k = 0; k < size; ++ k )
            {
                A[ j ][ k ] -= ( alpha * A[ index1 ][ k ] + gamma * A[ index2 ][ k ] );
                std::print( "j,k=({},{}) ", j, k );
                std::print( "index1,k=({},{}) ", index1, k );
                std::print( "index2,k=({},{}) \n", index2, k );
            }
            std::println();
            F[ j ] -= ( alpha * F[ index1 ] + gamma * F[ index2 ] );
        }
    }

    Print( A, "matrix A after cyclic reduction" );
    Print( F, "vector F after cyclic reduction" );
    //Back substitution
    int index = ( size - 1 ) / 2;
    x[ index ] = F[ index ] / A[ index ][ index ];
    std::print( "F[{}] = {}, A[{}][{}]={}\n", index, F[ index ], index, index, A[ index ][ index ] );
    std::print( "x[{}] = {}\n", index, x[ index ] );
    for ( int i = std::log2( size + 1 ) - 2; i >= 0; -- i )
    {
        std::print( "Back substitution i = {}\n", i );
        for ( int j = std::pow( 2, i + 1 ) - 1; j < size; j = j + std::pow( 2, i + 1 ) )
        {
            int offset = std::pow( 2, i );
            int index1 = j - offset;
            int index2 = j + offset;
            std::print( "index1,j,index2 = ({},{},{})\n", index1, j, index2 );
            x[ index1 ] = F[ index1 ];
            x[ index2 ] = F[ index2 ];
            for ( int k = 0; k < size; ++ k )
            {
                if ( k != index1 )
                {
                    x[ index1 ] -= A[ index1 ][ k ] * x[ k ];
                }
                if ( k != index2 )
                {
                    x[ index2 ] -= A[ index2 ][ k ] * x[ k ];
                }
            }
            x[ index1 ] = x[ index1 ] / A[ index1 ][ index1 ];
            x[ index2 ] = x[ index2 ] / A[ index2 ][ index2 ];
            std::print( "x[{}] = {}, x[{}] = {}\n", index1, x[ index1 ], index2, x[ index2 ] );
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

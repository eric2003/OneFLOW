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

void Print000( std::vector<std::vector<double>> & a, const std::string & name = "matrix" )
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

void Print( std::vector<std::vector<double>> & a, const std::string & name = "matrix" )
{
    int NI = a.size();
    std::cout << name << " = \n";
    for ( int i = 0; i < NI; ++ i )
    {
        int NJ = a[ i ].size();
        for ( int j = 0; j < NJ; ++ j )
        {
            std::cout << std::format( "{:5.2f} ", a[ i ][ j ] );
        }
        std::cout << "\n";
    }
    std::cout << "\n";
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
    //const int N = 3;
    const int N = 7;
    //const int N = 15;
    std::vector<double> a( N );
    std::vector<double> b( N );
    std::vector<double> c( N );
    std::vector<double> y( N );
    std::vector<double> x( N );
    std::vector<double> xx( N );
    for ( int i = 0; i < N; ++ i )
    {
        xx[ i ] = 0.0;
    }
    std::vector<double> F( N );
    std::vector<std::vector<double>> A( N );
    for ( int i = 0; i < N; ++ i )
    {
        A[ i ].resize( N );
        for ( int j = 0; j < N; ++ j )
        {
            A[ i ][ j ] = 0.0;
        }
        a[ i ] = 1;
        b[ i ] = -2;
        c[ i ] = 1;
        y[ i ] = i;
        F[ i ] = i;
    }
    A[ 0 ][ 0 ] = -2.0;
    A[ 0 ][ 1 ] = 1.0;
    A[ N - 1 ][ N - 2 ] = 1.0;
    A[ N - 1 ][ N - 1 ] = -2.0;

    for ( int i = 1; i < N - 1; ++ i )
    {
        A[ i ][ i ] = -2.0;
        A[ i ][ i - 1 ] = 1.0;
        A[ i ][ i + 1 ] = 1.0;
    }

    Print( A, "matrix A 0" );
    Print( F, "vector F" );

    std::vector<std::vector<double>> AA = A;
    std::vector<double> FF = F;

    //Cyclic reduction
    for ( int level = 0; level < std::log2( N + 1 ) - 1; ++ level )
    {
        std::print( "Cyclic reduction level = {}\n", level );
        int joffset = std::pow( 2, level );
        int jstep = std::pow( 2, level + 1 );
        int jstart = jstep - 1;
        for ( int j = jstart; j < N; j = j + jstep )
        {
            int in = j - joffset;
            int ip = j + joffset;
            double alpha = A[ j ][ in ] / A[ in ][ in ];
            double gamma = A[ j ][ ip ] / A[ ip ][ ip ];

            double alpha1 = a[ j ] / b[ in ];
            double gamma1 = c[ j ] / b[ ip ];
            std::print( "in,j,ip = ({},{},{})\n", in, j, ip );
            std::print( "A[{}][{}]={}\n", j, in, A[ j ][ in ] );
            std::print( "A[{}][{}]={}\n", j, ip, A[ j ][ ip ] );
            std::print( "A[{}][{}]={}\n", in, in, A[ in ][ in ] );
            std::print( "A[{}][{}]={}\n", ip, ip, A[ ip ][ ip ] );
            std::print( "alpha = {}, gamma = {}\n", alpha, gamma );
            for ( int k = 0; k < N; ++ k )
            {
                A[ j ][ k ] -= ( alpha * A[ in ][ k ] + gamma * A[ ip ][ k ] );
                std::print( "j,k=({},{}) ", j, k );
                std::print( "in,k=({},{}) ", in, k );
                std::print( "ip,k=({},{}) \n", ip, k );
            }
            std::println();
            b[ j ] -= ( alpha1 * c[ in ] + gamma1 * a[ ip ] );
            a[ j ] = - alpha1 * a[ in ];
            c[ j ] = - gamma1 * c[ ip ];
            y[ j ] -= ( alpha1 * y[ in ] + gamma1 * y[ ip ] );
            F[ j ] -= ( alpha * F[ in ] + gamma * F[ ip ] );
        }
    }

    Print( A, "matrix A after cyclic reduction" );
    Print( F, "vector F after cyclic reduction" );
    //Back substitution
    int index = ( N - 1 ) / 2;
    x[ index ] = F[ index ] / A[ index ][ index ];
    xx[ index ] = y[ index ] / b[ index ];
    std::print( "F[{}] = {}, A[{}][{}]={}\n", index, F[ index ], index, index, A[ index ][ index ] );
    std::print( "x[{}] = {}\n", index, x[ index ] );
    for ( int level = std::log2( N + 1 ) - 2; level >= 0; -- level )
    {
        std::print( "Back substitution level = {}\n", level );
        int joffset = std::pow( 2, level );
        int jstep = std::pow( 2, level + 1 );
        int jstart = jstep - 1;
        for ( int j = jstart; j < N; j = j + jstep )
        {
            int in = j - joffset;
            int ip = j + joffset;
            std::print( "in,j,ip = ({},{},{})\n", in, j, ip );
            x[ in ] = F[ in ];
            x[ ip ] = F[ ip ];
            for ( int k = 0; k < N; ++ k )
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

            if ( in - joffset < 0 )
            {
                xx[ in ] = ( y[ in ] - c[ in ] * xx[ in + joffset ] ) /  b[ in ];
            }
            else
            {
                xx[ in ] = ( y[ in ] - a[ in ] * xx[ in - joffset ] - c[ in ] * xx[ in + joffset ] ) / b[ in ];
            }

            if ( ip + joffset >= N )
            {
                xx[ ip ] = ( y[ ip ] - a[ ip ] * xx[ ip - joffset ] ) / b[ ip ];
            }
            else
            {
                xx[ ip ] = ( y[ ip ] - a[ ip ] * xx[ ip - joffset ] - c[ ip ] * xx[ ip + joffset ] ) / b[ ip ];
            }
        }
    }

    Print( x, "x" );
    Print( xx, "xx" );

    std::vector<double> ff( F.size() );
    MatrixMultiply( AA, x, ff );
    Print( ff, "ff" );
    Print( FF, "FF" );

    return 0;
}

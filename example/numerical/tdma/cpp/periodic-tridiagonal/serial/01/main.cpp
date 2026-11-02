import std;

template<typename T>
void Print( std::vector<T> & x, const std::string & name = "vector" )
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
    int NI = a.size();
    for ( int i = 0; i < NI; ++ i )
    {
        int NJ = a[ i ].size();
        for ( int j = 0; j < NJ; ++ j )
        {
            std::print( "{:4.1f} ", a[ i ][ j ] );
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

void FillA( const std::vector<double> &a, 
    const std::vector<double> &b, 
    const std::vector<double> &c, 
    std::vector<std::vector<double>> &A )
{
    int N = a.size();
    A.resize( N );
    for ( int i = 0; i < N; ++ i )
    {
        A[ i ].resize( N );
        if ( i >= 1 )
        {
            A[ i ][ i - 1 ] = a[ i ];
        }
        A[ i ][ i ] = b[ i ];
        if ( i < N - 1 )
        {
            A[ i ][ i + 1 ] = c[ i ];
        }
    }
}

void thomas_algorithm( const std::vector<double> & a,
    const std::vector<double> & b,
    const std::vector<double> & c,
    const std::vector<double> & d,
    std::vector<double> & x )
{
    size_t N = d.size();

    std::vector<double> c_star( N, 0.0 );
    std::vector<double> d_star( N, 0.0 );

    c_star[ 0 ] = c[ 0 ] / b[ 0 ];
    d_star[ 0 ] = d[ 0 ] / b[ 0 ];

    for ( int i = 1; i < N; ++ i )
    {
        double coef = 1.0 / ( b[ i ] - a[ i ] * c_star[ i - 1 ] );
        c_star[ i ] = c[ i ] * coef;
        d_star[ i ] = ( d[ i ] - a[ i ] * d_star[ i - 1 ] ) * coef;
    }

    x[ N - 1 ] = d_star[ N - 1 ];

    for ( int i = N - 2; i >= 0; -- i )
    {
        x[ i ] = d_star[ i ] - c_star[ i ] * x[ i + 1 ];
    }
}

void cyclic_thomas_algorithm( const std::vector<double> & a,
    const std::vector<double> & b,
    const std::vector<double> & c,
    const std::vector<double> & d,
    std::vector<double> & x )
{
    size_t N = a.size();
    std::vector<double> bNew = b;
    double gamma = - b[ 0 ];
    double aa = a[ 0 ] / gamma;
    bNew[ 0 ] -= gamma;
    bNew[ N - 1 ] -= c[ N - 1 ] * aa;
    std::vector<double> u( N, 0 );
    std::vector<double> v( N, 0 );
    u[ 0 ] = gamma;
    u[ N - 1 ] = c[ N - 1 ];

    v[ 0 ] = 1;
    v[ N - 1 ] = a[ 0 ] / gamma;

    std::vector<double> y( N, 0 );
    std::vector<double> q( N, 0 );

    thomas_algorithm( a, bNew, c, d, y );
    thomas_algorithm( a, bNew, c, u, q );
    double t1 = y[ 0 ] + aa * y[ N - 1 ];
    double t2 = q[ 0 ] + aa * q[ N - 1 ] + 1;
    double tt = t1 / t2;

    for ( int i = 0; i < N; ++ i )
    {
        x[ i ] = y[ i ] - tt * q[ i ];
    }
}
                                                                                                                                                                                             
int main( int argc, char ** argv )
{
    std::vector<std::vector<double>>A{ {2,1,0,7},{-1,7,4,0},{0,2,-3,2},{6,0,1,8} };
    Print( A, "A" );

    std::vector<double> q{ 32, 25, -5, 41 };

    std::vector<double> a{ 7, -1, 2, 1 };
    std::vector<double> b{ 2, 7, -3, 8 };
    std::vector<double> c{ 1, 4, 2,  6 };
    std::vector<double> d{ 32, 25, -5, 41 };

    std::vector<double> x( q.size() );

    cyclic_thomas_algorithm( a, b, c, d, x );
    Print( x, "x" );

    std::vector<double> y( d.size() );

    MatrixMultiply( A, x, y );
    Print( y, "y" );
    Print( d, "d" );

    return 0;
}

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

void ThomasAlgorithm( int N, double * b, double * a, double * c, double * x, double * q )
{
    double *l,*u,*d,*y;
    l = new double[ N ];
    u = new double[ N ];
    d = new double[ N ];
    y = new double[ N ];
    /* LU Decomposition */
    d[ 0 ] = a[ 0 ];
    u[ 0 ] = c[ 0 ];
    for ( int i = 1; i < N - 1; i++ )
    {
        l[ i ] = b[ i ] / d[ i - 1 ];
        d[ i ] = a[ i ] - l[ i ] * u[ i - 1 ];
        u[ i ] = c[ i ];
    }
    l[ N - 1 ] = b[ N - 1 ] / d[ N - 2 ];
    d[ N - 1 ] = a[ N - 1 ] - l[ N - 1 ] * u[ N - 2 ];
    /* Forward Substitution [L][y] = [q] */
    y[ 0 ] = q[ 0 ];
    for ( int i = 1; i < N; i++ )
    {
        y[ i ] = q[ i ] - l[ i ] * y[ i - 1 ];
    }
    /* Backward Substitution [U][x] = [y] */
    x[ N - 1 ] = y[ N - 1 ] / d[ N - 1 ];
    for ( int i = N - 2; i >= 0; i-- )
    {
        x[ i ] = ( y[ i ] - u[ i ] * x[ i + 1 ] ) / d[ i ];
    }
    delete[] l;
    delete[] u;
    delete[] d;
    delete[] y;
    return;
}

void ThomasAlgorithmLU(int N, double *b, double *a, double *c,
    double *l, double *u, double *d)
{
    /* LU Decomposition */
    d[ 0 ] = a[ 0 ];
    u[ 0 ] = c[ 0 ];
    for ( int i = 1; i < N -1; i++ )
    {
        l[ i ] = b[ i ] / d[ i - 1 ];
        d[ i ] = a[ i ] - l[ i ] * u[ i - 1 ];
        u[ i ] = c[ i ];
    }
    l[ N - 1 ] = b[ N - 1 ] / d[ N - 2 ];
    d[ N - 1 ] = a[ N - 1 ] - l[ N - 1 ] * u[ N - 2 ];
    return;
}

void ThomasAlgorithmSolve(int N, double *l, double *u, double *d,
    double *x, double *q)
{
    double * y = new double[ N ];
    /* Forward Substitution [L][y] = [q] */
    y[ 0 ] = q[ 0 ];
    for ( int i = 1; i < N; i++ )
    {
        y[ i ] = q[ i ] - l[ i ] * y[ i - 1 ];
    }
    /* Backward Substitution [U][x] = [y] */
    x[ N - 1 ] = y[ N - 1 ] / d[ N - 1 ];
    for ( int i = N - 2; i >= 0; i-- )
    {
        x[ i ] = ( y[ i ] - u[ i ] * x[ i + 1 ] ) / d[ i ];
    }
    delete[] y;
    return;
}

void thomas_algorithm( const std::vector<double> & a,
    const std::vector<double> & b,
    const std::vector<double> & c,
    std::vector<double> & d )
{
    size_t N = d.size();

    std::vector<double> c_star( N, 0.0 );
    std::vector<double> d_star( N, 0.0 );

    c_star[ 0 ] = c[ 0 ] / b[ 0 ];
    d_star[ 0 ] = d[ 0 ] / b[ 0 ];

    for ( int i = 1; i < N; ++ i )
    {
        double r = 1.0 / ( b[ i ] - a[ i ] * c_star[ i - 1 ] );
        d_star[ i ] = r * ( d[ i ] - a[ i ] * d_star[ i - 1 ] );
        c_star[ i ] = r * c[ i ];
    }

    d[ N - 1 ] = d_star[ N - 1 ];

    for ( int i = N - 2; i >= 0; -- i )
    {
        d[ i ] = d_star[ i ] - c_star[ i ] * d[ i + 1 ];
    }
}
                                                                                                                                                                                             
int main( int argc, char ** argv )
{
    std::vector<double> a{ 0, -1, -1, -1, -1 };
    std::vector<double> b{ 2,  2,  2,  2,  2 };
    std::vector<double> c{ -1, -1, -1, -1, 0 };
    std::vector<double> q{ 1.0, 1.0, 1.0, 1.0, 1.0 };

    std::vector<double> x( q.size() );
    std::vector<std::vector<double>> AA;
    FillA( a, b, c, AA );
    Print( AA, "AA" );

    int N = a.size();
    ThomasAlgorithm( N, a.data(),b.data(), c.data(), x.data(), q.data() );
    Print( x, "x" );

    std::vector<double> y( q.size() );
    MatrixMultiply( AA, x, y );
    Print( y, "y" );

    std::vector<double> l( q.size() );
    std::vector<double> u( q.size() );
    std::vector<double> d( q.size() );
    std::vector<double> xx( q.size() );

    ThomasAlgorithmLU( N, a.data(), b.data(), c.data(), l.data(), u.data(), d.data() );
    ThomasAlgorithmSolve( N, l.data(), u.data(), d.data(), xx.data(), q.data() );
    Print( xx, "xx" );
        
    return 0;
}

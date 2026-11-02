import std;

void Print( std::vector<double> & x, const std::string &name="vector")
{
    std::print( "{} = ", name );
    for ( auto v: x )
    {
        std::print( "{} ", v );
    }
    std::println();
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

void Boundary( std::vector<double> & a,
    std::vector<double> & b,
    std::vector<double> & c,
    std::vector<double> & y )
{
    int totalN = a.size();
    int N = totalN - 2;
    a[ 0 ] = 0;
    a[ 1 ] = 0;
    b[ 0 ] = 1;
    b[ N + 1 ] = 1;
    c[ N ] = 0;
    c[ N + 1 ] = 0;
    y[ 0 ] = 0;
    y[ N + 1 ] = 0;
}

void SetCyclicReductionValue
(
    std::vector<double> & a,
    std::vector<double> & b,
    std::vector<double> & c,
    std::vector<double> & y,
    std::vector<double> & a1,
    std::vector<double> & b1,
    std::vector<double> & c1,
    std::vector<double> & y1,
    std::vector<double> & a2,
    std::vector<double> & b2,
    std::vector<double> & c2,
    std::vector<double> & y2
)
{
    int N = a.size();
    int halfN = N / 2;
    a1.resize( halfN );
    b1.resize( halfN );
    c1.resize( halfN );
    y1.resize( halfN );
    a2.resize( halfN );
    b2.resize( halfN );
    c2.resize( halfN );
    y2.resize( halfN );
    for ( int i = 0; i < halfN; ++ i )
    {
        int ieven = 2 * i;
        int iodd  = ieven + 1;
        a1[ i ] = a[ ieven ];
        b1[ i ] = b[ ieven ];
        c1[ i ] = c[ ieven ];
        y1[ i ] = y[ ieven ];

        a2[ i ] = a[ iodd ];
        b2[ i ] = b[ iodd ];
        c2[ i ] = c[ iodd ];
        y2[ i ] = y[ iodd ];
    }
}

void CyclicReduction( std::vector<double> & a,
    std::vector<double> & b,
    std::vector<double> & c,
    std::vector<double> & y,
    std::vector<double> & abar,
    std::vector<double> & bbar,
    std::vector<double> & cbar,
    std::vector<double> & ybar
    )
{
    int totalN = a.size();
    int N = totalN - 2;
    for ( int i = 1; i <= N; ++ i )
    {
        int in = std::max( 0, i - 1 );
        int ip = std::min( N + 1, i + 1 );

        double alpha = - a[ i ] / b[ in ];
        double beta = - c[ i ] / b[ ip ];

        abar[ i ] = alpha * a[ in ];
        cbar[ i ] = beta * c[ ip ];
        bbar[ i ] = b[ i ] + alpha * c[ in ] + beta * a[ ip ];
        ybar[ i ] = y[ i ] + alpha * y[ in ] + beta * y[ ip ];
    }
}

void CR
(
    std::vector<double> & a,
    std::vector<double> & b,
    std::vector<double> & c,
    std::vector<double> & y,
    std::vector<double> & a1,
    std::vector<double> & b1,
    std::vector<double> & c1,
    std::vector<double> & y1,
    std::vector<double> & a2,
    std::vector<double> & b2,
    std::vector<double> & c2,
    std::vector<double> & y2
)
{
    int totalN = a.size();
    int N = totalN - 2;
    int halfN = N / 2;
    int totalHalfN = halfN + 2;
    a1.resize( totalHalfN );
    b1.resize( totalHalfN );
    c1.resize( totalHalfN );
    y1.resize( totalHalfN );

    a2.resize( totalHalfN );
    b2.resize( totalHalfN );
    c2.resize( totalHalfN );
    y2.resize( totalHalfN );
    for ( int i = 1; i <= halfN; ++ i )
    {
        int iodd = 2 * i - 1;
        int ieven = iodd + 1;

        a1[ i ] = a[ iodd ];
        b1[ i ] = b[ iodd ];
        c1[ i ] = c[ iodd ];
        y1[ i ] = y[ iodd ];

        a2[ i ] = a[ ieven ];
        b2[ i ] = b[ ieven ];
        c2[ i ] = c[ ieven ];
        y2[ i ] = y[ ieven ];
    }
}

void Create( std::vector<double> & a, std::vector<double> & aNew )
{
    aNew.assign( a.begin() + 1, a.end() - 1 );
}

void crtridiag( std::vector<double> & a,
    std::vector<double> & b,
    std::vector<double> & c,
    std::vector<double> & y,
    std::vector<double> & x
    )
{
    int totalN = a.size();
    int N = totalN - 2;
    x.resize( N );
    if ( N == 1 )
    {
        x[ 0 ] = y[ 1 ] / b[ 1 ];
        return;
    }

    Boundary( a, b, c, y );

    std::vector<double> acr( a.size() );
    std::vector<double> bcr( b.size() );
    std::vector<double> ccr( c.size() );
    std::vector<double> ycr( y.size() );

    CyclicReduction( a, b, c, y, acr, bcr, ccr, ycr );

    std::vector<double> a1, b1, c1, y1;
    std::vector<double> a2, b2, c2, y2;

    CR( acr, bcr, ccr, ycr, a1, b1, c1, y1, a2, b2, c2, y2 );

    std::vector<double> x1;
    crtridiag( a1, b1, c1, y1, x1 );

    Print( x1, std::format( "x1(N={})", N ) );
    std::vector<double> x2;
    crtridiag( a2, b2, c2, y2, x2 );
    Print( x2, std::format( "x2(N={})", N ) );

    std::vector<double> xx;
    for ( int i = 0; i < x1.size(); ++ i )
    {
        xx.push_back( x1[ i ] );
        xx.push_back( x2[ i ] );
    }
    Print( xx, std::format( "x1+2(N={})", N ) );
    x = xx;
}
                                                                                                                                                                                             
int main( int argc, char ** argv )
{
    const int N = 8;

    int totalN = N + 2;
    std::vector<double> a( totalN, -1 );
    std::vector<double> b( totalN, 2 );
    std::vector<double> c( totalN, -1 );
    std::vector<double> y( totalN, 1.0 );
    std::vector<double> x;

    crtridiag( a, b, c, y, x );

    Print( x, "xFinal" );

    return 0;
}

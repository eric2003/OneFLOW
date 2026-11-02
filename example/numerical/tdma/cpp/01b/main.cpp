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

void split(
    std::vector<double> & a,
    std::vector<double> & b,
    std::vector<double> & c,
    std::vector<double> & y,
    std::vector<double> & a_bar,
    std::vector<double> & b_bar,
    std::vector<double> & c_bar,
    std::vector<double> & y_bar )
{
    int N = a.size();
    a_bar.resize( N );
    b_bar.resize( N );
    c_bar.resize( N );
    y_bar.resize( N );
    for ( int i = 0; i < N; ++ i )
    {
        double bim1 = 0;
        double bip1 = 0;
        double aim1 = 0;
        double aip1 = 0;
        double cim1 = 0;
        double cip1 = 0;
        double yim1 = 0;
        double yip1 = 0;
        double alpha = 0;
        double beta = 0;
        if ( i != 0 )
        {
            aim1 =  a[ i - 1 ];
            bim1 =  b[ i - 1 ];
            cim1 =  c[ i - 1 ];
            yim1 =  y[ i - 1 ];
            alpha = - a[ i ] / bim1;
        }

        if ( i != N - 1 )
        {
            aip1 =  a[ i + 1 ];
            bip1 =  b[ i + 1 ];
            cip1 =  c[ i + 1 ];
            yip1 =  y[ i + 1 ];
            beta = - c[ i ] / bip1;
        }
        a_bar[ i ] = alpha * aim1;
        c_bar[ i ] = beta * cip1;
        b_bar[ i ] = b[ i ] + alpha * cim1 + beta * aip1;
        y_bar[ i ] = y[ i ] + alpha * yim1 + beta * yip1;
    }
}
                                                                                                                                                                                             
int main( int argc, char ** argv )
{
    const int N = 4;

    std::vector<double> a( N, -1 );
    std::vector<double> b( N, 2 );
    std::vector<double> c( N, -1 );
    std::vector<double> d( N, 0.0 );
    a[ 0 ] = 0;
    c[ N - 1 ] = 0;
    d[ 0 ] = 1;
    d[ N - 1 ] = 1;

    std::vector<double> x( d.size() );

    std::vector<double> as, bs, cs, ds;
    split( a, b, c, d, as, bs, cs, ds );

    Print( a, "a" );
    Print( b, "b" );
    Print( c, "c" );
    Print( d, "d" );

    Print( as, "as" );
    Print( bs, "bs" );
    Print( cs, "cs" );
    Print( ds, "ds" );

    //thomas_algorithm( a, b, c, d, x );

    //for ( auto v : x )
    //{
    //    std::print( "{} ", v );
    //}
 
    return 0;
}

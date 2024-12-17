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
                                                                                                                                                                                             
int main( int argc, char ** argv )
{
    int p = 3;
    int N0 = std::pow( 2, p ) - 1;
    int N = N0 + 5;
    int totalN = std::pow( 2, p + 1 ) - 1;
    std::vector<double> a( totalN, -1 );
    std::vector<double> b( totalN, 2 );
    std::vector<double> c( totalN, -1 );
    std::vector<double> y( totalN, 0.0 );
    a[ 0 ] = 0;
    c[ N - 1 ] = 0;
     
    y[ 0 ] = 1;
    y[ N - 1 ] = 1;

    for ( int i = N; i < totalN; ++ i )
    {
        a[ i ] = 0;
        c[ i ] = 0;
        b[ i ] = 1;
        y[ i ] = 0;
    }
     
    //for ( int i = 0; i < N; ++ i )
    //{
    //    y[ i ] = i;
    //}
     
    std::vector<double> x( y.size() );
     
    thomas_algorithm( a, b, c, y, x );
    Print( x, "x" );

    int kkk = 1;

    return 0;
}

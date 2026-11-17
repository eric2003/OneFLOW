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
    for ( int N = 4; N >= 1; --N )
    {
        int totalN = N + 2;
        std::vector<double> a( totalN, -1 );
        std::vector<double> b( totalN, 2 );
        std::vector<double> c( totalN, -1 );
        std::vector<double> y( totalN, 0.0 );
        a[ 0 ] = 0;
        a[ 1 ] = 0;
        b[ 0 ] = 1;
        b[ N + 1 ] = 1;
        c[ N ] = 0;
        c[ N + 1 ] = 0;
        y[ 0 ] = 0;
        y[ 1 ] = 1;
        y[ N ] = 1;
        y[ N + 1 ] = 0;

        std::vector<double> anew( a.begin() + 1, a.end() - 1 );
        std::vector<double> bnew( b.begin() + 1, b.end() - 1 );
        std::vector<double> cnew( c.begin() + 1, c.end() - 1 );
        std::vector<double> ynew( y.begin() + 1, y.end() - 1 );
        std::vector<double> x( ynew.size() );

        thomas_algorithm( anew, bnew, cnew, ynew, x );
        std::print( "N={}\n", N );
        Print( x, "x" );
    }
    {
        int N = 8;
        double lv = std::log2( 8 );
        int M = log2( N + 1 ) - 1;
        std::print( "M={}\n", M );
        int offset = 1;
        int joffset = 2;
        int start = 0;
        for ( int i = 0; i <= M; ++ i )
        {
            int jstart = joffset - 1;
            std::print( "i={}\n", i );
            std::print( "jstart={},N-1={}\n", jstart, N - 1 );

            std::print( "j=" );

            for ( int j = jstart; j <= N - 1; j += joffset )
            {
                int id1 = j - offset;
                int id2 = j + offset;
                std::print( "{} ", j );
            }
            std::print( "\n" );

            std::print( "id1=" );

            for ( int j = jstart; j <= N - 1; j += joffset )
            {
                int id1 = j - offset;
                int id2 = j + offset;
                std::print( "{} ", id1 );
            }
            std::print( "\n" );

            std::print( "id2=" );

            for ( int j = jstart; j <= N - 1; j += joffset )
            {
                int id1 = j - offset;
                int id2 = j + offset;
                std::print( "{} ", id2 );
            }
            std::print( "\n" );

            offset *= 2;
            joffset *= 2;
        }
        int kkk = 1;

    }

    return 0;
}

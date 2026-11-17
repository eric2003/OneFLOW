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
    //const int N = 3;
    const int N = 7;
    //const int N = 15;
    std::vector<double> a( N );
    std::vector<double> b( N );
    std::vector<double> c( N );
    std::vector<double> y( N );
    std::vector<double> x( N );
    for ( int i = 0; i < N; ++ i )
    {
        a[ i ] = 1;
        b[ i ] = -2;
        c[ i ] = 1;
        y[ i ] = i;
    }

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

            double alpha1 = a[ j ] / b[ in ];
            double gamma1 = c[ j ] / b[ ip ];
            std::print( "in,j,ip = ({},{},{})\n", in, j, ip );
            std::println();
            b[ j ] -= ( alpha1 * c[ in ] + gamma1 * a[ ip ] );
            a[ j ] = - alpha1 * a[ in ];
            c[ j ] = - gamma1 * c[ ip ];
            y[ j ] -= ( alpha1 * y[ in ] + gamma1 * y[ ip ] );
        }
    }

    //Back substitution
    int index = ( N - 1 ) / 2;
    x[ index ] = y[ index ] / b[ index ];
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

            if ( in - joffset < 0 )
            {
                x[ in ] = ( y[ in ] - c[ in ] * x[ in + joffset ] ) /  b[ in ];
            }
            else
            {
                x[ in ] = ( y[ in ] - a[ in ] * x[ in - joffset ] - c[ in ] * x[ in + joffset ] ) / b[ in ];
            }

            if ( ip + joffset >= N )
            {
                x[ ip ] = ( y[ ip ] - a[ ip ] * x[ ip - joffset ] ) / b[ ip ];
            }
            else
            {
                x[ ip ] = ( y[ ip ] - a[ ip ] * x[ ip - joffset ] - c[ ip ] * x[ ip + joffset ] ) / b[ ip ];
            }
        }
    }

    Print( x, "x" );

    return 0;
}

import std;

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
        double r = 1.0 / ( b[ i ] - a[ i ] * c_star[ i - 1 ] );
        d_star[ i ] = r * ( d[ i ] - a[ i ] * d_star[ i - 1 ] );
        c_star[ i ] = r * c[ i ];
    }

    x[ N - 1 ] = d_star[ N - 1 ];

    for ( int i = N - 2; i >= 0; -- i )
    {
        x[ i ] = d_star[ i ] - c_star[ i ] * x[ i + 1 ];
    }
}
                                                                                                                                                                                             
int main( int argc, char ** argv )
{
    std::vector<double> a{ 0, -1, -1, -1, -1 }; // 下对角线  
    std::vector<double> b{ 2,  2,  2,  2,  2 }; // 主对角线  
    std::vector<double> c{ -1, -1, -1, -1, 0 }; // 上对角线  
    std::vector<double> d{ 1.0, 1.0, 1.0, 1.0, 1.0 }; // 右边的常数向量  
    std::vector<double> x( d.size() ); // 结果向量  
                                                                                                                                                                      
    thomas_algorithm( a, b, c, d, x );
    for ( auto v : x )
    {
        std::print( "{} ", v );
    }
 
    return 0;
}

#include <vector>  
#include <print>

void Print( std::vector<double> & x );
void Print( std::vector<std::vector<double>> & a );
void MatrixMultiply( std::vector<std::vector<double>> & a, std::vector<double> & x, std::vector<double> & y );
void gaussElimination( std::vector<std::vector<double>> & a, std::vector<double> & b );

void Print( std::vector<double> & x )
{
    for ( auto v: x )
    {
        std::print( "{} ", v );
    }
    std::println();
}

void Print( std::vector<std::vector<double>> & a )
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

void gaussElimination( std::vector<std::vector<double>> & a, std::vector<double> & b )
{
    int n = a.size();  

    // 消元过程  
    for ( int i = 0; i < n; ++ i ) {
        // 找到当前列的最大值  
        double maxVal = std::abs( a[ i ][ i ] );
        int maxRow = i;  
        for ( int k = i + 1; k < n; ++ k )
        {
            if ( std::abs( a[ k ][ i ] ) > maxVal )
            {
                maxVal = std::abs( a[ k ][ i ] );
                maxRow = k;
            }  
        }  
        std::swap( a[ i ], a[ maxRow ] );
        std::swap( b[ i ], b[ maxRow ] );

        // 进行消元  
        for ( int j = i + 1; j < n; ++ j )
        {
            double factor = a[ j ][ i ] / a[ i ][ i ];
            for ( int k = i; k < n; ++ k )
            {
                a[ j ][ k ] -= factor * a[ i ][ k ];
            }  
            b[ j ] -= factor * b[ i ];
        }  
    }  

    // 回代过程  
    std::vector<double> x( n );
    for ( int i = n - 1; i >= 0; -- i )
    {
        x[ i ] = b[ i ];
        for ( int j = i + 1; j < n; ++ j )
        {
            x[ i ] -= a[ i ][ j ] * x[ j ];
        }  
        x[ i ] /= a[ i ][ i ];
    }  

    // 输出结果  
    for ( int i = 0; i < n; ++ i )
    {  
        std::print( "b{} = {:12.6f}\n", i + 1, b[ i ] );
    }  

    // 输出结果  
    for ( int i = 0; i < n; ++ i )
    {  
        std::print( "x{} = {:12.6f}\n", i + 1, x[ i ] );
    }  
}  

int main() 
{  
    std::vector<std::vector<double>>a{
        {3.0, 2.0,-4.0},
        {2.0, 3.0, 3.0},
        {5.0, -3, 1.0}
    };

    Print( a );

    //std::print( "a.size()={}\n", a.size() );
    std::vector<double> b{ 3.0,15.0,14.0 };
    std::vector<double> xx{ 3.0,1.0,2.0 };
    std::vector<double> y( xx.size() );

    MatrixMultiply( a, xx, y );

    gaussElimination( a, b );

    Print( a );

    Print( y );

    return 0;
} 


import std;

void upper_triangular_matrix( const std::vector<std::vector<double>> & a )
{
}
                                                                                                                                                                                             
int main( int argc, char ** argv )
{
    std::vector<std::vector<double>>u{
        {10,11,12},
        {0,20,21},
        {0,0,30}
    };
    int N = u.size();
    for ( int i = 0; i < N; ++ i )
    {
        for ( int j = 0; j < N; ++ j )
        {
            std::print( "{} ", u[ i ][ j ] );
        }
        std::println();
    }

    return 0;
}

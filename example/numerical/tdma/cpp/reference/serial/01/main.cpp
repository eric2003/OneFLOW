import std;

const int size = 15;
                                                                                                                                                                                   
int main( int argc, char ** argv )
{
    int j,k;
    int index1, index2, offset;
    double alpha, gamma;
    /* Part 1 */
    double * x = new double[ size ];
    for ( int i = 0; i < size; ++ i )
    {
        x[ i ] = 0.0;
    }
    double * F = new double[ size ];
    double ** A = new double * [ size ];
    for ( int i = 0; i < size; ++ i )
    {
        A[ i ] = new double[ size ];
        for ( int j = 0; j < size; ++ j )
        {
            A[ i ][ j ] = 0.0;
        }
        F[ i ] = (double)i;
    }
    A[ 0 ][ 0 ] = -2.0;
    A[ 0 ][ 1 ] = 1.0;
    A[ size - 1 ][ size - 2 ] = 1.0;
    A[ size - 1 ][ size - 1 ] = -2.0;
    for ( int i = 1; i < size - 1; ++ i )
    {
        A[ i ][ i ] = -2.0;
        A[ i ][ i - 1 ] = 1.0;
        A[ i ][ i + 1 ] = 1.0;
    }

    delete [] x;
    delete [] F;
    for ( int i = 0; i < size; ++ i )
    {
        delete[] A[ i ];
    }
    delete [] A;

    return 0;
}

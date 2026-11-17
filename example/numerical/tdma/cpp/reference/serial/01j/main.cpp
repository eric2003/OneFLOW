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

void Print000( std::vector<std::vector<double>> & a, const std::string & name = "matrix" )
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

void Print( std::vector<std::vector<double>> & a, const std::string & name = "matrix" )
{
    int NI = a.size();
    std::cout << name << " = \n";
    for ( int i = 0; i < NI; ++ i )
    {
        int NJ = a[ i ].size();
        for ( int j = 0; j < NJ; ++ j )
        {
            std::cout << std::format( "{:5.2f} ", a[ i ][ j ] );
        }
        std::cout << "\n";
    }
    std::cout << "\n";
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

void SetProcessMatrixValue( std::vector<std::vector<std::vector<double>>> &AList,int N, int nProc )
{
    AList.resize( nProc );
    for ( int iProc = 0; iProc < nProc; ++ iProc )
    {
        int nRow = 3;
        AList[ iProc ].resize( nRow );
        for ( int i = 0; i < nRow; ++ i )
        {
            AList[ iProc ][ i ].resize( N, 0 );
        }

    }
    for ( int iProc = 0; iProc < nProc; ++ iProc )
    {
        std::vector<std::vector<double>> & A = AList[ iProc ];
        //int loc = 2 * iProc + 1; //1,3,5,7,9,11,13
        int loc = 2 * iProc; //0,2,4,6,8,10,12
        std::vector<int> global_locs( 3 );
        for ( int i = 0; i < 3; ++ i )
        {
            int iloc = loc + i;
            global_locs[ i ] = iloc;
            if ( ( iloc - 1 ) >= 0 )
            {
                A[ i ][ iloc - 1 ] =  1.0;
            }
            A[ i ][ iloc     ] = -2.0;
            if ( ( iloc + 1 ) < N )
            {
                A[ i ][ iloc + 1 ] =  1.0;
            }
        }
        std::print("iProc={}\n", iProc );
        Print( global_locs, "global_locs" );
        Print( A, "matrix A" );
    }
}


int main( int argc, char ** argv )
{
   //p =   0, N =     0, nprocs =    -1
   //p =   1, N =     1, nprocs =     0
   //p =   2, N =     3, nprocs =     1
   //p =   3, N =     7, nprocs =     3
   //p =   4, N =    15, nprocs =     7
   //p =   5, N =    31, nprocs =    15
   //p =   6, N =    63, nprocs =    31
   //p =   7, N =   127, nprocs =    63
   //p =   8, N =   255, nprocs =   127
   //p =   9, N =   511, nprocs =   255
   //p =  10, N =  1023, nprocs =   511
    //Memory allocation and generation of matrix

    int nprocs = 3;
    int N = std::pow( 2, std::log2( nprocs + 1 ) + 1 ) - 1;

    std::vector<double> a( N );
    std::vector<double> b( N );
    std::vector<double> c( N );
    std::vector<double> y( N );
    std::vector<double> x( N );
    std::vector<double> xx( N );
    for ( int i = 0; i < N; ++ i )
    {
        xx[ i ] = 0.0;
    }
    std::vector<double> F( N );
    std::vector<std::vector<double>> A( N );
    for ( int i = 0; i < N; ++ i )
    {
        A[ i ].resize( N );
        for ( int j = 0; j < N; ++ j )
        {
            A[ i ][ j ] = 0.0;
        }
        a[ i ] = 1;
        b[ i ] = -2;
        c[ i ] = 1;
        y[ i ] = i;
        F[ i ] = i;
    }
    A[ 0 ][ 0 ] = -2.0;
    A[ 0 ][ 1 ] = 1.0;
    A[ N - 1 ][ N - 2 ] = 1.0;
    A[ N - 1 ][ N - 1 ] = -2.0;

    for ( int i = 1; i < N - 1; ++ i )
    {
        A[ i ][ i ] = -2.0;
        A[ i ][ i - 1 ] = 1.0;
        A[ i ][ i + 1 ] = 1.0;
    }

    Print( A, "matrix A 0" );
    Print( F, "vector F" );

    std::vector<std::vector<std::vector<double>>> AProcs;
    SetProcessMatrixValue( AProcs, N, nprocs );

    std::vector<std::vector<double>> AA = A;
    std::vector<double> FF = F;

    //Cyclic reduction
    for ( int level = 0; level < std::log2( N + 1 ) - 1; ++ level )
    {
        std::print( "Cyclic reduction level = {}\n", level );
        int joffset = std::pow( 2, level );
        int jstep = std::pow( 2, level + 1 );
        int jstart = jstep - 1;
        std::print( "j = " );
        for ( int j = jstart; j < N; j = j + jstep )
        {
            std::print( "{} ", j );
        }
        std::println();
        std::print( "in = " );
        for ( int j = jstart; j < N; j = j + jstep )
        {
            int in = j - joffset;
            std::print( "{} ", in );
        }
        std::println();
        std::print( "ip = " );
        for ( int j = jstart; j < N; j = j + jstep )
        {
            int ip = j + joffset;
            std::print( "{} ", ip );
        }
        std::println();

        //j=1,3,5
        for ( int j = jstart; j < N; j = j + jstep )
        {
            int in = j - joffset;
            int ip = j + joffset;
            double alpha = A[ j ][ in ] / A[ in ][ in ];
            double gamma = A[ j ][ ip ] / A[ ip ][ ip ];

            double alpha1 = a[ j ] / b[ in ];
            double gamma1 = c[ j ] / b[ ip ];
            std::print( "in,j,ip = ({},{},{})\n", in, j, ip );
            std::print( "alpha = {}, gamma = {}\n", alpha, gamma );
            for ( int k = 0; k < N; ++ k )
            {
                A[ j ][ k ] -= ( alpha * A[ in ][ k ] + gamma * A[ ip ][ k ] );
            }
            std::println();
            b[ j ] -= ( alpha1 * c[ in ] + gamma1 * a[ ip ] );
            a[ j ] = - alpha1 * a[ in ];
            c[ j ] = - gamma1 * c[ ip ];
            y[ j ] -= ( alpha1 * y[ in ] + gamma1 * y[ ip ] );
            F[ j ] -= ( alpha * F[ in ] + gamma * F[ ip ] );
        }

        std::vector<int> jlist;
        jlist.push_back(1);
        jlist.push_back(3);
        jlist.push_back(5);
        std::vector<int> jlocals;
        jlocals.push_back(1);
        jlocals.push_back(1);
        jlocals.push_back(1);

        for ( int iProc = 0; iProc < nprocs; ++ iProc )
        {
            std::vector<std::vector<double>> &AProc = AProcs[iProc];
            Print( AProc, std::format( "iProc = {} Matrix AProc 000 ", iProc ) );
        }

        for ( int iProc = 0; iProc < nprocs; ++ iProc )
        {
            std::vector<std::vector<double>> &AProc = AProcs[iProc];
            int j = jlist[iProc];
            int jl = jlocals[iProc];

            int inl = jl - joffset;
            int ipl = jl + joffset;
            int in = j - joffset;
            int ip = j + joffset;
            double alpha = AProc[ jl ][ in ] / AProc[ inl ][ in ];
            double gamma = AProc[ jl ][ ip ] / AProc[ ipl ][ ip ];	
            for ( int k = 0; k < N; ++ k )
            {
                AProc[ jl ][ k ] -= ( alpha * AProc[ inl ][ k ] + gamma * AProc[ ipl ][ k ] );
            }
            Print( AProc, std::format( "level = {} iProc = {} Matrix AProc ", level, iProc ) );
        }
        int kkk = 1;
    }

    Print( A, "matrix A after cyclic reduction" );
    Print( F, "vector F after cyclic reduction" );
    //Back substitution
    int index = ( N - 1 ) / 2;
    x[ index ] = F[ index ] / A[ index ][ index ];
    xx[ index ] = y[ index ] / b[ index ];
    std::print( "F[{}] = {}, A[{}][{}]={}\n", index, F[ index ], index, index, A[ index ][ index ] );
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
            x[ in ] = F[ in ];
            x[ ip ] = F[ ip ];
            for ( int k = 0; k < N; ++ k )
            {
                if ( k != in )
                {
                    x[ in ] -= A[ in ][ k ] * x[ k ];
                }
                if ( k != ip )
                {
                    x[ ip ] -= A[ ip ][ k ] * x[ k ];
                }
            }
            x[ in ] = x[ in ] / A[ in ][ in ];
            x[ ip ] = x[ ip ] / A[ ip ][ ip ];
            std::print( "x[{}] = {}, x[{}] = {}\n", in, x[ in ], ip, x[ ip ] );

            if ( in - joffset < 0 )
            {
                xx[ in ] = ( y[ in ] - c[ in ] * xx[ in + joffset ] ) /  b[ in ];
            }
            else
            {
                xx[ in ] = ( y[ in ] - a[ in ] * xx[ in - joffset ] - c[ in ] * xx[ in + joffset ] ) / b[ in ];
            }

            if ( ip + joffset >= N )
            {
                xx[ ip ] = ( y[ ip ] - a[ ip ] * xx[ ip - joffset ] ) / b[ ip ];
            }
            else
            {
                xx[ ip ] = ( y[ ip ] - a[ ip ] * xx[ ip - joffset ] - c[ ip ] * xx[ ip + joffset ] ) / b[ ip ];
            }
        }
    }

    Print( x, "x" );
    Print( xx, "xx" );

    std::vector<double> ff( F.size() );
    MatrixMultiply( AA, x, ff );
    Print( ff, "ff" );
    Print( FF, "FF" );

    return 0;
}

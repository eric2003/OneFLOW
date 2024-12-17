import std;

void Print( std::vector<std::vector<double>> & a, const std::string & name = "matrix" )
{
	int N = a.size();
	for ( int i = 0; i < N; ++ i )
	{
		for ( int j = 0; j < N; ++ j )
		{
			std::print( "{:2.0f} ", a[ i ][ j ] );
		}
		std::println();
	}
	std::println();
}

int main( int argc, char ** argv )
{
	int N = 15;
	int NProc = 7;
	std::vector<std::vector<double>> A;
	A.resize( N );
	for ( int i = 0; i < N; ++ i )
	{
		A[ i ].resize( N );
	}
	for ( int iProc = 0; iProc < NProc; ++ iProc )
	{
		//int loc = 2 * iProc + 1; //1,3,5,7,9,11,13
		int loc = 2 * iProc; //0,2,4,6,8,10,12
		for ( int i = 0; i < 3; ++ i )
		{
			int iloc = loc + i;
			if ( ( iloc - 1 ) >= 0 )
			{
				A[ iloc ][ iloc - 1 ] =  1.0;
			}
			A[ iloc ][ iloc     ] = -2.0;
			if ( ( iloc + 1 ) < N )
			{
				A[ iloc ][ iloc + 1 ] =  1.0;
			}
		}
	}

	Print( A, "matrix A" );
    return 0;
}

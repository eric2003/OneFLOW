import std;

void Print( std::vector<std::vector<double>> & a, const std::string & name = "matrix" )
{
	int NI = a.size();
	for ( int i = 0; i < NI; ++ i )
	{
		int NJ = a[ i ].size();
		for ( int j = 0; j < NJ; ++ j )
		{
			std::print( "{:2.0f} ", a[ i ][ j ] );
		}
		std::println();
	}
	std::println();
}

void SetMatrixValue( int N, int nProc )
{
	std::vector<std::vector<double>> A;
	A.resize( N );
	for ( int i = 0; i < N; ++ i )
	{
		A[ i ].resize( N );
	}
	for ( int iProc = 0; iProc < nProc; ++ iProc )
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
}

void SetMatrixValueNew( int N, int nProc )
{
	std::vector<std::vector<std::vector<double>>> AList;
	AList.resize( nProc );
	for ( int iProc = 0; iProc < nProc; ++ iProc )
	{
		int nRow = 3;
		AList[ iProc ].resize( nRow );
		for ( int i = 0; i < nRow; ++ i )
		{
			AList[ iProc ][ i ].resize( N );
		}

	}
	for ( int iProc = 0; iProc < nProc; ++ iProc )
	{
		std::vector<std::vector<double>> & A = AList[ iProc ];
		//int loc = 2 * iProc + 1; //1,3,5,7,9,11,13
		int loc = 2 * iProc; //0,2,4,6,8,10,12
		for ( int i = 0; i < 3; ++ i )
		{
			int iloc = loc + i;
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
		Print( A, "matrix A" );
	}
}

int main( int argc, char ** argv )
{
	int N = 15;
	int nProc = 7;

	SetMatrixValue( N, nProc );
	SetMatrixValueNew( N, nProc );

    return 0;
}

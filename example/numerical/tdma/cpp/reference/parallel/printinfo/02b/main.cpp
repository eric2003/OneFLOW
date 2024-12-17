import std;

int main( int argc, char ** argv )
{
	for ( int p = 0; p <= 10; ++ p )
	{
		int N = std::pow( 2, p ) - 1;
		int nprocs = ( N + 1 ) / 2 - 1;
		std::print( "p = {:3}, N = {:5}, nprocs = {:5}\n", p, N, nprocs );
	}

    return 0;
}

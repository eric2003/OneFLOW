import std;

int main( int argc, char ** argv )
{
	for ( int nprocs = 1; nprocs <= 7; ++ nprocs )
	{
		int N = (int)std::pow( 2, std::log2( nprocs + 1 ) + 1 ) - 1;
		std::print( "nprocs = {}, N = {:4}, (N-1)/2 = {}\n", nprocs, N, ( N - 1 ) / 2 );
	}

    return 0;
}

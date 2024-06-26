#include <stdio.h>
#ifdef ENABLE_OPENACC
#include <openacc.h>
#endif

#include "myopenacc.h"
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "laplace2d.h"


#ifdef ENABLE_OPENACC
void test_open_acc()
{
    const int N = 5;
    int a[N];
    #pragma acc parallel loop
    for ( int i = 0; i < N; ++ i )
    {
        a[i] = i;
    }
    
    #pragma acc parallel
    for ( int i = 0; i < N; ++ i )
    {
        printf( "openacc a[%d] = %d\n", i, a[i] );
    }
    std::cout << "Hello from myopenacc.cpp using std::cout\n";
	
}
#endif


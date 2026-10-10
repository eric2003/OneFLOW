#include <stdio.h>
#include "jacobi.h"
#include <iostream>
#include "laplace2d.h"
#include <string.h>
#include <stdio.h>
#include <vector>
#include <omp.h>

void Jacobi_Test()
{
    const int n = 4096;
    const int m = 4096;
    const int iter_max = 1000;
    
    const double tol = 1.0e-6;
    double error = 1.0;

    // RAII: vector owns the buffers; pass .data() to the C-style kernels.
    std::vector< double > A( static_cast< size_t >( n ) * static_cast< size_t >( m ) );
    std::vector< double > Anew( static_cast< size_t >( n ) * static_cast< size_t >( m ) );
    
    initialize( A.data(), Anew.data(), m, n );
        
    printf("Jacobi relaxation Calculation: %d x %d mesh\n", n, m);
    
    double st = omp_get_wtime();
    int iter = 0;
   
    while ( error > tol && iter < iter_max )
    {
        error = calcNext( A.data(), Anew.data(), m, n );
        swap( A.data(), Anew.data(), m, n );

        if(iter % 100 == 0) printf("%5d, %0.6f\n", iter, error);
        
        iter++;

    }

    double runtime = omp_get_wtime() - st;
 
    printf(" total: %f s\n", runtime);

    // A / Anew destroyed automatically; do not free() vector storage.
}


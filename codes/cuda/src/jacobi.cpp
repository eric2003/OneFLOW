#include <stdio.h>
#include "jacobi.h"
#include <iostream>
#include "laplace2d.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
using namespace std;

void Jacobi_Test()
{
    const int n = 4096;
    const int m = 4096;
    const int iter_max = 1000;
    
    const double tol = 1.0e-6;
    double error = 1.0;

    double *restrict A    = (double*)malloc(sizeof(double)*n*m);
    double *restrict Anew = (double*)malloc(sizeof(double)*n*m);
    
    initialize(A, Anew, m, n);
        
    printf("Jacobi relaxation Calculation: %d x %d mesh\n", n, m);
    
    double st = omp_get_wtime();
    int iter = 0;
   
    while ( error > tol && iter < iter_max )
    {
        error = calcNext(A, Anew, m, n);
        swap(A, Anew, m, n);

        if(iter % 100 == 0) printf("%5d, %0.6f\n", iter, error);
        
        iter++;

    }

    double runtime = omp_get_wtime() - st;
 
    printf(" total: %f s\n", runtime);

    deallocate(A, Anew);

}


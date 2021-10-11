#ifdef ENABLE_CUDA

#include <stdio.h>
#include <cuda_runtime.h>
#include "cuda_sub.h"
#include <iostream>
using namespace std;

void GetCudaDeviceCount( int &num_gpus )
{
    cudaGetDeviceCount( &num_gpus );

    for ( int i = 0; i < num_gpus; ++ i )
    {
        cudaDeviceProp dprop;
        cudaGetDeviceProperties(&dprop, i);
        printf("%d: %s\n", i, dprop.name);
    }
}

#endif

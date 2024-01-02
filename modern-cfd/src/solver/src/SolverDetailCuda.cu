/*---------------------------------------------------------------------------*\
OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
Copyright (C) 2017-2024 He Xin and the OneFLOW contributors.
-------------------------------------------------------------------------------
License
This file is part of OneFLOW.

OneFLOW is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

OneFLOW is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
for more details.

You should have received a copy of the GNU General Public License
along with OneFLOW.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#ifdef PRJ_ENABLE_CUDA
#include "Cmpi.h"
#include <omp.h>
#include <cstdio>
#include <cuda_runtime.h>

void SolverInitCuda()
{
    cudaGetDeviceCount( &Cmpi::num_gpus );
    if ( Cmpi::num_gpus < 1 ) {
        std::printf("no CUDA capable devices were detected\n");
        //std::exit(1);
    }

    std::printf("number of host CPUs:\t%d\n", omp_get_num_procs());
    std::printf("number of CUDA devices:\t%d\n", Cmpi::num_gpus);

    for ( int i = 0; i < Cmpi::num_gpus; ++ i )
    {
        cudaDeviceProp dprop;
        cudaGetDeviceProperties( &dprop, i);
        std::printf("   %d: %s\n", i, dprop.name);
    }

    std::printf("---------------------------\n");

    int nCpuThreads = 8;
    omp_set_num_threads( nCpuThreads );
#ifdef ENABLE_OPENMP
#pragma omp parallel
#endif
    {
        unsigned int cpu_thread_id = omp_get_thread_num();
        unsigned int num_cpu_threads = omp_get_num_threads();
        //std::printf( "Solver::Solver() CPU thread %d (of %d)\n", cpu_thread_id, num_cpu_threads );
    }
}

void SetDeviceCuda( int cpu_thread_id )
{
    int gpu_id = -1;
    cudaSetDevice( cpu_thread_id % Cmpi::num_gpus );
    cudaGetDevice( &gpu_id );
}

__global__ void GpuCfdCopyVector( float *a, const float *b, int numElements )
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if ( i < numElements )
    {
        a[i] = b[i];
    }
}

__global__ void GpuCfdScalarUpdate( float * q, const float * qn, float c, const float * timestep, const float * ds, int ni )
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if ( i < ni + 1 && i > 0 )
    {
        float cfl = c * timestep[ i ] / ds[ i ];
        q[ i ] = qn[ i ] - cfl * ( qn[ i ] - qn[ i - 1 ] );
    }
}

void CfdCopyVectorCuda( float * a, float * b, int ni )
{
    std::size_t nSize = ni * sizeof(float);

    float * dev_a;
    float * dev_b;
    cudaMalloc((void **)&dev_a, nSize);
    cudaMalloc((void **)&dev_b, nSize);

    cudaMemcpy(dev_a, a, nSize, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, b, nSize, cudaMemcpyHostToDevice);

    int block_size = 256;
    int grid_size = ( ni + block_size - 1 ) / block_size;
    dim3 grid_dim( grid_size );
    dim3 block_dim( block_size );  // 256 threads per block

    GpuCfdCopyVector<<<grid_dim, block_dim>>>( dev_a, dev_b, ni );
    cudaDeviceSynchronize();
    cudaMemcpy(a, dev_a, nSize, cudaMemcpyDeviceToHost);
    cudaFree(dev_a);
    cudaFree(dev_b);
}

void CfdScalarUpdateCuda( float * q, float * qn, float c, float * timestep, float * ds, int ni )
{
    float * dev_q;
    float * dev_qn;
    float * dev_timestep;
    float * dev_ds;
    int nElem = ni + 2;
    std::size_t nSize = nElem * sizeof(float);

    cudaMalloc((void **)&dev_qn, nSize);
    cudaMalloc((void **)&dev_q, nSize);
    cudaMalloc((void **)&dev_timestep, nSize);
    cudaMalloc((void **)&dev_ds, nSize);

    cudaMemcpy(dev_qn, qn, nSize, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_q, q, nSize, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_ds, ds, nSize, cudaMemcpyHostToDevice);
    cudaMemcpy(dev_timestep, timestep, nSize, cudaMemcpyHostToDevice);

    int block_size = 256;
    int grid_size = ( nElem + block_size - 1 ) / block_size;
    dim3 grid_dim( grid_size );
    dim3 block_dim( block_size );  // 256 threads per block

    //std::printf("Solver::SolveField CUDA kernel launch with %d blocks of %d threads\n", grid_size, block_size);
    GpuCfdScalarUpdate<<<grid_dim, block_dim>>>(dev_q, dev_qn, c, dev_timestep, dev_ds, ni);
    cudaDeviceSynchronize();
    cudaMemcpy(q, dev_q, nSize, cudaMemcpyDeviceToHost);
    cudaFree(dev_q);
    cudaFree(dev_qn);
    cudaFree(dev_timestep);
    cudaFree(dev_ds);
}
#endif

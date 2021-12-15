/*---------------------------------------------------------------------------*\
    OneFLOW - LargeScale Multiphysics Scientific Simulation Environment
    Copyright (C) 2017-2021 He Xin and the OneFLOW contributors.
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

#include "AAASolverCuda.h"
#include "HXMath.h"
#include "Constant.h"
#ifdef ENABLE_CUDA
#include <cuda_runtime.h>
#endif
#include <iostream>

BeginNameSpace( ONEFLOW )

__global__ void SetValueKernel(Real *dev_a, Real *dev_b, int *dev_id );
__global__ void addKernel(int *a, int *b, int *c );

__global__ void SetValueKernel(Real *dev_a, Real *dev_b, int *dev_id )
{
    int iface = threadIdx.x;
    int icell = dev_id[ iface ];
    dev_a[ iface ] = dev_b[ icell ];
}

__global__ void SetValueKernelReal(Real *dev_a, Real *dev_b, int *dev_id, int nFaces, int nTCells )
{
    int iface = blockDim.x * blockIdx.x + threadIdx.x;
    if ( iface < nFaces )
    {
        int jcell = dev_id[ iface ];
        dev_a[iface] = dev_b[jcell];
    }
}

__global__ void MyInvFluxCuda(Real * qf1, Real * qf2, Real * invflux, Real * xfn, Real * yfn, Real * zfn, Real * area, int nFaces )
{
    int iFace = blockDim.x * blockIdx.x + threadIdx.x;
    if ( iFace < nFaces )
    {
        Real vxl = 1.0;
        Real vyl = 0.0;
        Real vzl = 0.0;

        Real vxr = 1.0;
        Real vyr = 0.0;
        Real vzr = 0.0;

        Real q_L = qf1[ iFace ];
        Real q_R = qf2[ iFace ];

        Real vnl  = xfn[ iFace ] * vxl + yfn[ iFace ] * vyl + zfn[ iFace ] * vzl;
        Real vnr  = xfn[ iFace ] * vxr + yfn[ iFace ] * vyr + zfn[ iFace ] * vzr;

        Real eigenL = vnl;
        Real eigenR = vnr;

        //eigenL = half * ( eigenL + ABS( eigenL ) );
        //eigenR = half * ( eigenR - ABS( eigenR ) );
        eigenL = 0.5 * ( eigenL + abs( eigenL ) );
        eigenR = 0.5 * ( eigenR - abs( eigenR ) );

        Real fL = q_L * eigenL;
        Real fR = q_R * eigenR;
        Real fM = fL + fR;

        Real areaM = area[ iFace ];
        invflux[ iFace ] = fM * areaM;
    }
}

__global__ void MyAddF2CFieldCudaDevice(Real * fField, Real * cField, int * lc, int * rc, int nBFaces, int nFaces )
{
    int iFace = blockDim.x * blockIdx.x + threadIdx.x;
    if ( iFace < nBFaces )
    {
        int lc_ = lc[ iFace ];
        Real value = fField[ iFace ];
        //cField[ lc_ ] -= fField[ iFace ];
        atomicAdd( &cField[ lc_ ], - value );
    }
    else if ( ( iFace >= nBFaces ) && ( iFace < nFaces ) )
    {
        int lc_ = lc[ iFace ];
        int rc_ = rc[ iFace ];

        Real value = fField[ iFace ];
        //cField[ lc_ ] -= fField[ iFace ];
        //cField[ rc_ ] += fField[ iFace ];
        atomicAdd( &cField[ lc_ ], - value );
        atomicAdd( &cField[ rc_ ], value );
    }
}

__global__ void MyAddF2CFieldCudaDeviceNoAtomic(Real * fField, Real * cField, int * lc, int * rc, int nBFaces, int nFaces )
{
    int iFace = blockDim.x * blockIdx.x + threadIdx.x;
    if ( iFace < nBFaces )
    {
        int lc_ = lc[ iFace ];
        cField[ lc_ ] -= fField[ iFace ];
    }
    else if ( ( iFace >= nBFaces ) && ( iFace < nFaces ) )
    {
        int lc_ = lc[ iFace ];
        int rc_ = rc[ iFace ];

        cField[ lc_ ] -= fField[ iFace ];
        cField[ rc_ ] += fField[ iFace ];
    }
}

__global__ void MyZoneTimeIntergralCudaDevice(Real * res, Real * vol, Real dt, int nCells)
{
    int iCell = blockDim.x * blockIdx.x + threadIdx.x;
    if ( iCell < nCells )
    {
        Real ovol = 1.0 / vol[ iCell ];
        Real coef = dt * ovol;
        res[ iCell ] *= coef;
    }
}

__global__ void MyZoneUpdateCudaDevice(Real *q, Real *res, int nCells)
{
    int iCell = blockDim.x * blockIdx.x + threadIdx.x;
    if ( iCell < nCells )
    {
        q[ iCell ] += res[ iCell ];
    }
}

void addWithCuda(int *a, int *b, int *c, unsigned int nElems)
{
    int *dev_a = 0;
    int *dev_b = 0;
    int *dev_c = 0;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaSetDevice(0);

    cudaMalloc((void**)&dev_a, nElems * sizeof(int));
    cudaMalloc((void**)&dev_b, nElems * sizeof(int));
    cudaMalloc((void**)&dev_c, nElems * sizeof(int));

    cudaMemcpy(dev_a, a, nElems * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, b, nElems * sizeof(int), cudaMemcpyHostToDevice);

    // Launch a kernel on the GPU with one thread for each element.
    addKernel<<<1, nElems>>>(dev_a, dev_b, dev_c);

    cudaDeviceSynchronize();

    cudaMemcpy(c, dev_c, nElems * sizeof(int), cudaMemcpyDeviceToHost);

    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);
}

__global__ void addKernel(int *a, int *b, int *c )
{
    int i = threadIdx.x;
    c[i] = a[i] + b[i];
}

__global__ void addRealKernel(Real *a, Real *b, Real *c )
{
    int i = threadIdx.x;
    c[i] = a[i] + b[i];
}

__global__ void addRealSwapKernel(Real *a, Real *b, int *id, Real *c )
{
    int i = threadIdx.x;
    int j = id[ i ];
    c[i] = a[i] + b[j];
}

__global__ void setRealSwapKernel(Real *a, int *id, Real *c )
{
    int i = threadIdx.x;
    int j = id[ i ];
    c[i] = a[j];
}

__global__ void setRealSwapKernelNew(Real *a, Real *b, int *id  )
{
    int i = threadIdx.x;
    int j = id[ i ];
    a[i] = b[j];
}

__global__ void setRealSwapKernelNew1(Real *a, Real *b, int *id, int nElems )
{
    //int i = threadIdx.x;
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if ( i < nElems )
    {
        int j = id[ i ];
        a[i] = b[j];
    }
}

void addRealWithCuda(Real *a, Real *b, Real *c, unsigned int nElems)
{
    Real *dev_a = 0;
    Real *dev_b = 0;
    Real *dev_c = 0;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaSetDevice(0);

    cudaMalloc((void**)&dev_a, nElems * sizeof(Real));
    cudaMalloc((void**)&dev_b, nElems * sizeof(Real));
    cudaMalloc((void**)&dev_c, nElems * sizeof(Real));

    cudaMemcpy(dev_a, a, nElems * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, b, nElems * sizeof(Real), cudaMemcpyHostToDevice);

    // Launch a kernel on the GPU with one thread for each element.
    addRealKernel<<<1, nElems>>>(dev_a, dev_b, dev_c);

    cudaDeviceSynchronize();

    cudaMemcpy(c, dev_c, nElems * sizeof(Real), cudaMemcpyDeviceToHost);

    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);
}

void addRealSwapWithCuda(Real *a, Real *b, int * id, Real *c, unsigned int nElems)
{
    Real *dev_a = 0;
    Real *dev_b = 0;
    Real *dev_c = 0;
    int * dev_id = 0;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaSetDevice(0);

    cudaMalloc((void**)&dev_a, nElems * sizeof(Real));
    cudaMalloc((void**)&dev_b, nElems * sizeof(Real));
    cudaMalloc((void**)&dev_c, nElems * sizeof(Real));
    cudaMalloc((void**)&dev_id, nElems * sizeof(Real));

    cudaMemcpy(dev_a, a, nElems * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, b, nElems * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_id, id, nElems * sizeof(int), cudaMemcpyHostToDevice);

    // Launch a kernel on the GPU with one thread for each element.
    addRealSwapKernel<<<1, nElems>>>(dev_a, dev_b, dev_id, dev_c);

    cudaDeviceSynchronize();

    cudaMemcpy(c, dev_c, nElems * sizeof(Real), cudaMemcpyDeviceToHost);

    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_id);
}

void setRealSwapWithCuda(Real *a, int * id, Real *c, unsigned int nElems)
{
    Real *dev_a = 0;
    Real *dev_c = 0;
    int * dev_id = 0;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaSetDevice(0);

    cudaMalloc((void**)&dev_a, nElems * sizeof(Real));
    cudaMalloc((void**)&dev_c, nElems * sizeof(Real));
    cudaMalloc((void**)&dev_id, nElems * sizeof(Real));

    cudaMemcpy(dev_a, a, nElems * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_id, id, nElems * sizeof(int), cudaMemcpyHostToDevice);

    // Launch a kernel on the GPU with one thread for each element.
    setRealSwapKernel<<<1, nElems>>>(dev_a, dev_id, dev_c);

    cudaDeviceSynchronize();

    cudaMemcpy(c, dev_c, nElems * sizeof(Real), cudaMemcpyDeviceToHost);

    cudaFree(dev_c);
    cudaFree(dev_a);
    cudaFree(dev_id);
}

void setRealSwapWithCudaNew(Real *a, Real *b, int * id,  unsigned int nElems)
{
    Real *dev_a = 0;
    Real *dev_b = 0;
    int * dev_id = 0;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaSetDevice(0);

    cudaMalloc((void**)&dev_a, nElems * sizeof(Real));
    cudaMalloc((void**)&dev_b, nElems * sizeof(Real));
    cudaMalloc((void**)&dev_id, nElems * sizeof(Real));

    cudaMemcpy(dev_b, b, nElems * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_id, id, nElems * sizeof(int), cudaMemcpyHostToDevice);

    //setRealSwapKernelNew<<<1, nElems>>>(dev_a, dev_b, dev_id);

    int threadsPerBlock = 256;
    int blocksPerGrid =(nElems + threadsPerBlock - 1) / threadsPerBlock;

    printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);

    setRealSwapKernelNew1<<<blocksPerGrid, threadsPerBlock>>>(dev_a, dev_b, dev_id, nElems);

    cudaDeviceSynchronize();

    cudaMemcpy(a, dev_a, nElems * sizeof(Real), cudaMemcpyDeviceToHost);

    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_id);
}

void setRealSwapWithCudaNewRealProblem(Real *a, Real *b, int * id, unsigned int nFaces, unsigned int nCells)
{
    Real *dev_a = 0;
    Real *dev_b = 0;
    int * dev_id = 0;

    // Choose which GPU to run on, change this on a multi-GPU system.
    cudaSetDevice(0);

    cudaMalloc((void**)&dev_a, nFaces * sizeof(Real));
    cudaMalloc((void**)&dev_b, nCells * sizeof(Real));
    cudaMalloc((void**)&dev_id, nFaces * sizeof(Real));

    cudaMemcpy(dev_b, b, nCells * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_id, id, nFaces * sizeof(int), cudaMemcpyHostToDevice);

    setRealSwapKernelNew<<<1, nFaces>>>(dev_a, dev_b, dev_id);

    cudaDeviceSynchronize();

    cudaMemcpy(a, dev_a, nFaces * sizeof(Real), cudaMemcpyDeviceToHost);

    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_id);
}

void SetValueWithCuda(Real *aface, Real *bcell, int *id, unsigned int nFaces, unsigned int nTCells)
{
    Real *dev_a = 0;
    Real *dev_b = 0;
    int *dev_id = 0;

    cudaSetDevice(0);

    cudaMalloc((void**)&dev_a, nFaces * sizeof(Real));
    cudaMalloc((void**)&dev_b, nTCells * sizeof(Real));
    cudaMalloc((void**)&dev_id, nFaces * sizeof(int));

    cudaMemcpy(dev_a, aface, nFaces * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, bcell, nTCells * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_id, id, nFaces * sizeof(int), cudaMemcpyHostToDevice);

    int threadsPerBlock = 256;
    int blocksPerGrid = (nFaces + threadsPerBlock - 1) / threadsPerBlock;

    //printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);

    SetValueKernelReal<<<blocksPerGrid, threadsPerBlock>>>(dev_a, dev_b, dev_id, nFaces, nTCells );

    cudaDeviceSynchronize();

    cudaMemcpy(aface, dev_a, nFaces * sizeof(Real), cudaMemcpyDeviceToHost);

    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_id);
}

void MyCalcInvFluxCuda(Real *qf1, Real *qf2, Real *invflux, Real *xfn, Real *yfn, Real *zfn, Real *area, int nFaces)
{
    Real *dev_qf1 = 0;
    Real *dev_qf2 = 0;
    Real *dev_invflux = 0;
    Real *dev_xfn = 0;
    Real *dev_yfn = 0;
    Real *dev_zfn = 0;
    Real *dev_area = 0;

    cudaSetDevice(0);

    cudaMalloc((void**)&dev_qf1, nFaces * sizeof(Real));
    cudaMalloc((void**)&dev_qf2, nFaces * sizeof(Real));
    cudaMalloc((void**)&dev_invflux, nFaces * sizeof(Real));
    cudaMalloc((void**)&dev_xfn, nFaces * sizeof(Real));
    cudaMalloc((void**)&dev_yfn, nFaces * sizeof(Real));
    cudaMalloc((void**)&dev_zfn, nFaces * sizeof(Real));
    cudaMalloc((void**)&dev_area, nFaces * sizeof(Real));

    cudaMemcpy(dev_qf1, qf1, nFaces * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_qf2, qf2, nFaces * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_invflux, invflux, nFaces * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_xfn, xfn, nFaces * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_yfn, yfn, nFaces * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_zfn, zfn, nFaces * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_area, area, nFaces * sizeof(Real), cudaMemcpyHostToDevice);


    int threadsPerBlock = 256;
    int blocksPerGrid = (nFaces + threadsPerBlock - 1) / threadsPerBlock;

    //printf("CUDA kernel launch with %d blocks of %d threads\n", blocksPerGrid, threadsPerBlock);

    MyInvFluxCuda<<<blocksPerGrid, threadsPerBlock>>>(dev_qf1, dev_qf2, dev_invflux, dev_xfn, dev_yfn, dev_zfn, dev_area, nFaces );

    cudaDeviceSynchronize();

    cudaMemcpy(invflux, dev_invflux, nFaces * sizeof(Real), cudaMemcpyDeviceToHost);

    cudaFree(dev_qf1);
    cudaFree(dev_qf2);
    cudaFree(dev_invflux);
    cudaFree(dev_xfn);
    cudaFree(dev_yfn);
    cudaFree(dev_zfn);
    cudaFree(dev_area);
}

void MyAddF2CFieldCuda(Real *fField, Real *cField, int *lc, int * rc, int nBFaces, int nFaces, int nTCells)
{
    Real *dev_fField = 0;
    Real *dev_cField = 0;
    int *dev_lc = 0;
    int *dev_rc = 0;

    cudaSetDevice(0);

    cudaMalloc((void**)&dev_fField, nFaces * sizeof(Real));
    cudaMalloc((void**)&dev_cField, nTCells * sizeof(Real));
    cudaMalloc((void**)&dev_lc, nFaces * sizeof(int));
    cudaMalloc((void**)&dev_rc, nFaces * sizeof(int));

    cudaMemcpy(dev_fField, fField, nFaces * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_cField, cField, nTCells * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_lc, lc, nFaces * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_rc, rc, nFaces * sizeof(int), cudaMemcpyHostToDevice);

    int threadsPerBlock = 256;
    int blocksPerGrid = (nFaces + threadsPerBlock - 1) / threadsPerBlock;

    MyAddF2CFieldCudaDevice<<<blocksPerGrid, threadsPerBlock>>>(dev_fField, dev_cField, dev_lc, dev_rc, nBFaces, nFaces );
    //MyAddF2CFieldCudaDeviceNoAtomic<<<blocksPerGrid, threadsPerBlock>>>(dev_fField, dev_cField, dev_lc, dev_rc, nBFaces, nFaces );

    cudaDeviceSynchronize();

    cudaMemcpy(cField, dev_cField, nTCells * sizeof(Real), cudaMemcpyDeviceToHost);

    cudaFree(dev_fField);
    cudaFree(dev_cField);
    cudaFree(dev_lc);
    cudaFree(dev_rc);
}


void MyZoneTimeIntergralCuda(Real *res, Real *vol, Real dt, int nCells)
{
    Real *dev_res = 0;
    Real *dev_vol = 0;

    cudaSetDevice(0);

    cudaMalloc((void**)&dev_res, nCells * sizeof(Real));
    cudaMalloc((void**)&dev_vol, nCells * sizeof(Real));

    cudaMemcpy(dev_res, res, nCells * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_vol, vol, nCells * sizeof(Real), cudaMemcpyHostToDevice);

    int threadsPerBlock = 256;
    int blocksPerGrid = (nCells + threadsPerBlock - 1) / threadsPerBlock;

    MyZoneTimeIntergralCudaDevice<<<blocksPerGrid, threadsPerBlock>>>(dev_res, dev_vol, dt, nCells);

    cudaDeviceSynchronize();

    cudaMemcpy(res, dev_res, nCells * sizeof(Real), cudaMemcpyDeviceToHost);

    cudaFree(dev_res);
    cudaFree(dev_vol);

}

void MyZoneUpdateCuda(Real *q, Real *res, int nCells)
{
    Real *dev_res = 0;
    Real *dev_q = 0;

    cudaSetDevice(0);

    cudaMalloc((void**)&dev_res, nCells * sizeof(Real));
    cudaMalloc((void**)&dev_q, nCells * sizeof(Real));

    cudaMemcpy(dev_res, res, nCells * sizeof(Real), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_q, q, nCells * sizeof(Real), cudaMemcpyHostToDevice);

    int threadsPerBlock = 256;
    int blocksPerGrid = (nCells + threadsPerBlock - 1) / threadsPerBlock;

    MyZoneUpdateCudaDevice<<<blocksPerGrid, threadsPerBlock>>>(dev_q, dev_res, nCells);

    cudaDeviceSynchronize();

    cudaMemcpy(q, dev_q, nCells * sizeof(Real), cudaMemcpyDeviceToHost);

    cudaFree(dev_q);
    cudaFree(dev_res);
}

EndNameSpace

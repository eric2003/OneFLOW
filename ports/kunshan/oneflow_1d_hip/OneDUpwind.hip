#include "OneDUpwind.h"
#include <hip/hip_runtime.h>
#include <cmath>
#include <stdexcept>
#include <string>

namespace {
__global__ void ReconstructKernel( const double *u, int start, int nx,
                                    double *left, double *right ) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if ( i > nx ) return;
    int ii = i - 1;
    left[i] = u[ii - start];
    right[i] = u[ii + 1 - start];
}
__global__ void ResidualKernel( int nx, double cp, double cn, double invdx,
                                const double *left, const double *right,
                                double *res ) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if ( i >= nx ) return;
    double fip = cp * left[i + 1] + cn * right[i + 1];
    double fim = cp * left[i] + cn * right[i];
    res[i] = -(fip - fim) * invdx;
}
void Check( hipError_t e, const char *name ) {
    if ( e != hipSuccess ) {
        throw std::runtime_error(
            std::string( name ) + ": " + hipGetErrorString( e ) );
    }
}
}

void OneDHipUpwindResidual( const double *u, int uSize, int uStart, int nx,
                            double c, double dx, double *left, double *right,
                            double *residual ) {
    double *du=nullptr,*dl=nullptr,*dr=nullptr,*ds=nullptr;
    Check( hipMalloc(&du, uSize*sizeof(double)), "hipMalloc u" );
    Check( hipMalloc(&dl, (nx+1)*sizeof(double)), "hipMalloc left" );
    Check( hipMalloc(&dr, (nx+1)*sizeof(double)), "hipMalloc right" );
    Check( hipMalloc(&ds, nx*sizeof(double)), "hipMalloc residual" );
    Check( hipMemcpy(du,u,uSize*sizeof(double),hipMemcpyHostToDevice), "copy u" );
    int faceBlocks = (nx + 1 + 255) / 256;
    int cellBlocks = (nx + 255) / 256;
    double cp=0.5*(c+std::abs(c)), cn=0.5*(c-std::abs(c));
    hipLaunchKernelGGL(ReconstructKernel, dim3(faceBlocks), dim3(256), 0, 0,
                       du,uStart,nx,dl,dr);
    Check( hipGetLastError(), "reconstruction kernel launch" );
    hipLaunchKernelGGL(ResidualKernel, dim3(cellBlocks), dim3(256), 0, 0,
                       nx,cp,cn,1.0/dx,dl,dr,ds);
    Check( hipGetLastError(), "residual kernel launch" );
    Check( hipDeviceSynchronize(), "kernel synchronize" );
    Check( hipMemcpy(left,dl,(nx+1)*sizeof(double),hipMemcpyDeviceToHost), "copy left" );
    Check( hipMemcpy(right,dr,(nx+1)*sizeof(double),hipMemcpyDeviceToHost), "copy right" );
    Check( hipMemcpy(residual,ds,nx*sizeof(double),hipMemcpyDeviceToHost), "copy residual" );
    hipFree(ds); hipFree(dr); hipFree(dl); hipFree(du);
}

#ifdef ENABLE_CUDA

#include <cuda_runtime.h>
#include <iostream>
using namespace std;

void printDeviceProp(const cudaDeviceProp& prop)
{
    std::cout << "Device Name : " << prop.name << "\n";
    std::cout << "totalGlobalMem : " << prop.totalGlobalMem << "\n";
    std::cout << "sharedMemPerBlock " << prop.sharedMemPerBlock << "\n";
    std::cout << "regsPerBlock : " << prop.regsPerBlock << "\n";
    std::cout << "warpSize :" << prop.warpSize << "\n";
    std::cout << "memPitch : " << prop.memPitch << "\n";
    std::cout << "maxThreadsPerBlock " << prop.maxThreadsPerBlock << "\n";
    std::cout << "maxThreadsDim[0 - 2] : " << prop.maxThreadsDim[0] << " " << prop.maxThreadsDim[1] << " " << prop.maxThreadsDim[2] << "\n";
    std::cout << "maxGridSize[0 - 2] " << prop.maxGridSize[0] << " " << prop.maxGridSize[1] << " " << prop.maxGridSize[2] << "\n";
    std::cout << "totalConstMem : " << prop.totalConstMem << "\n";
    std::cout << "major.minor : " << prop.major << "." << prop.minor << "\n";
    std::cout << "clockRate : " << prop.clockRate << "\n";
    std::cout << "textureAlignment :" << prop.textureAlignment << "\n";
    std::cout << "deviceOverlap : " << prop.deviceOverlap << "\n";
    std::cout << "multiProcessorCount : " << prop.multiProcessorCount << "\n";
}

bool InitCUDA()
{
    //used to count the device numbers
    int count;

    cudaGetDeviceCount(&count);
    if (count == 0) {
        fprintf(stderr, "There is no device.\n");
        return false;
    }

    // find the device >= 1.X
    int i;
    for (i = 0; i < count; ++i) {
        cudaDeviceProp prop;
        if (cudaGetDeviceProperties(&prop, i) == cudaSuccess) {
            if (prop.major >= 1) {
                printDeviceProp(prop);
                break;
            }
        }
    }

    // if can't find the device
    if (i == count) {
        fprintf(stderr, "There is no device supporting CUDA 1.x.\n");
        return false;
    }

    // set cuda device
    cudaSetDevice(i);

    return true;
}

#endif // ENABLE_CUDA

#pragma once

#ifdef ENABLE_CUDA

#include <cuda_runtime.h>
void printDeviceProp(const cudaDeviceProp& prop);
bool InitCUDA();

#endif

#pragma once

#ifdef ENABLE_CUDA

#include <cuda_runtime.h>
void printDeviceProp( const cudaDeviceProp & prop, int device_idx );
bool InitCUDA();

#endif

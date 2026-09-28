#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#ifdef ENABLE_CUDA

void GetCudaDeviceCount( int & num_gpus );

#endif

#ifdef __cplusplus
} // closing brace for extern "C"
#endif

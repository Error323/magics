#ifndef GAMMA_CORRECT_CUH
#define GAMMA_CORRECT_CUH

#include <cuda_runtime_api.h>

namespace gpu
{
__device__ float sinc(float x);
__device__ float L(float x);

__global__ void transpose(uchar4 *odata,
                          uchar4 *idata,
                          int width,
                          int height);

__global__ void gamma(uchar4 *src,
                      uchar4 *dst,
                      float g,
                      int n);

__global__ void lanczos(uchar4 *src,
                        uchar4 *dst,
                        int nSrc,
                        int nDst,
                        float factor,
                        float scale,
                        float support);

bool process(uchar4 *src,
             const int srcSize,
             uchar4 *dst,
             const int dstSize,
             const float g);
}

#endif // GAMMA_CORRECT_CUH

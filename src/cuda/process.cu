#include "process.h"
#include "init.h"
#include "../Verbose.hpp"
#include "../Debugger.hpp"
#include "../../crawler/Constants.hpp"

#include <math.h>
#include <math_constants.h>

#define BLOCK_DIM_2D 16
#define MAX_SCALE 64.0f
#define MAX_KERNEL_SIZE int(MAX_SCALE*LANCZOS_WIDTH*2.0f)

namespace gpu
{
int iDivUp(int a, int b)
{
  return (a % b != 0) ? (a / b + 1) : (a / b);
}

__device__ float sinc(float x)
{
  x *= CUDART_PI_F;
  return sinf(x) / x;
}

__device__ float L(float x)
{
  if (x == 0.0f)
    return 1.0f;

  return sinc(x) * sinc(x/LANCZOS_WIDTH);
}

__device__ float4 rgba2float(uchar4 p)
{
  float4 rgba;
  rgba.x = p.x / 255.0f;
  rgba.y = p.y / 255.0f;
  rgba.z = p.z / 255.0f;
  rgba.w = p.w / 255.0f;
  return rgba;
}

__device__ uchar4 pixel(float4 p)
{
  uchar4 rgba;
  rgba.x = min(max(p.x, 0.0f), 255.0f) + 0.5f;
  rgba.y = min(max(p.y, 0.0f), 255.0f) + 0.5f;
  rgba.z = min(max(p.z, 0.0f), 255.0f) + 0.5f;
  rgba.w = 255;
  return rgba;
}

__device__ uchar4 float2rgba(float4 p)
{
  uchar4 rgba;
  rgba.x = __saturatef(p.x) * 255.0f;
  rgba.y = __saturatef(p.y) * 255.0f;
  rgba.z = __saturatef(p.z) * 255.0f;
  rgba.w = __saturatef(p.w) * 255.0f;
  return rgba;
}



__global__ void gamma(uchar4 *src, uchar4 *dst, float g, int n)
{
  int x_i = blockIdx.x * BLOCK_DIM_2D + threadIdx.x;
  int y_i = blockIdx.y * BLOCK_DIM_2D + threadIdx.y;

  if (x_i >= n || y_i >= n)
    return;

  int index = y_i*n + x_i;
  float4 pixel = rgba2float(src[index]);
  pixel.x = powf(pixel.x, g);
  pixel.y = powf(pixel.y, g);
  pixel.z = powf(pixel.z, g);
  dst[index] = float2rgba(pixel);
}



__global__ void lanczos(uchar4 *src, uchar4 *dst, int nSrc, int nDst, float factor, float scale, float support)
{
  int x = blockIdx.x*BLOCK_DIM_2D + threadIdx.x;
  int y = blockIdx.y*BLOCK_DIM_2D + threadIdx.y;

  if (x >= nDst || y >= nDst)
    return;

  float center_y = (y + 0.5f) * scale;
  float center_x = (x + 0.5f) * scale;
  int start_y = max(int(center_y-support+0.5f), 0);
  int start_x = max(int(center_x-support+0.5f), 0);
  int stop_y  = min(int(center_y+support+0.5f), nSrc);
  int stop_x  = min(int(center_x+support+0.5f), nSrc);
  int n_y = stop_y - start_y;
  int n_x = stop_x - start_x;

  // compute kernels
  float kernel_x[MAX_KERNEL_SIZE];
  float kernel_y[MAX_KERNEL_SIZE];

  float density = 0.0f;
  for (int i = 0; i < n_y; i++)
  {
    float phase = start_y + i - center_y + 0.5f;
    density += kernel_y[i] = L(factor*phase);
  }
  for (int i = 0; i < n_y; i++)
    kernel_y[i] /= density;

  density = 0.0f;
  for (int i = 0; i < n_x; i++)
  {
    float phase = start_x + i - center_x + 0.5f;
    density += kernel_x[i] = L(factor*phase);
  }
  for (int i = 0; i < n_x; i++)
    kernel_x[i] /= density;

  // apply kernels and store result in resized image
  float4 p = make_float4(0.0f, 0.0f, 0.0f, 0.0f);
  uchar4 rgba;
  float lanczos_xy;
  for (int i = 0; i < n_y; i++)
  {
    for (int j = 0; j < n_x; j++)
    {
      lanczos_xy = kernel_y[i] * kernel_x[j];
      rgba = src[(i+start_y)*nSrc+(j+start_x)];
      p.x += lanczos_xy * rgba.x;
      p.y += lanczos_xy * rgba.y;
      p.z += lanczos_xy * rgba.z;
    }
  }
  dst[y*nDst+x] = pixel(p);
}



bool process(uchar4 *src, const int srcSize,
             uchar4 *dst, const int dstSize,
             const float g)
{
  size_t src_mem = srcSize*srcSize*sizeof(uchar4);
  size_t dst_mem = dstSize*dstSize*sizeof(uchar4);
  if (src_mem+dst_mem > gpu::freeMemory())
    return false;

  float factor  = dstSize / float(srcSize);
  float scale   = 1.0f / factor;
  float support = scale * LANCZOS_WIDTH;

  if (scale > MAX_SCALE)
    return false;

  Debug("[gpu]");
  dim3 blocks, threads;
  uchar4 *d_src, *d_dst;
  cudaSafeCall(cudaMalloc((void**)&d_src, src_mem));
  cudaSafeCall(cudaMalloc((void**)&d_dst, dst_mem));
  cudaSafeCall(cudaMemcpy(d_src, src, src_mem, cudaMemcpyHostToDevice));

  blocks  = dim3(iDivUp(srcSize, BLOCK_DIM_2D), iDivUp(srcSize, BLOCK_DIM_2D));
  threads = dim3(BLOCK_DIM_2D, BLOCK_DIM_2D);
  gamma<<<blocks, threads>>>(d_src, d_src, g, srcSize);

  blocks  = dim3(iDivUp(dstSize, BLOCK_DIM_2D), iDivUp(dstSize, BLOCK_DIM_2D));
  threads = dim3(BLOCK_DIM_2D, BLOCK_DIM_2D);
  lanczos<<<blocks, threads>>>(d_src, d_dst, srcSize, dstSize, factor, scale, support);

  blocks  = dim3(iDivUp(dstSize, BLOCK_DIM_2D), iDivUp(dstSize, BLOCK_DIM_2D));
  threads = dim3(BLOCK_DIM_2D, BLOCK_DIM_2D);
  gamma<<<blocks, threads>>>(d_dst, d_dst, 1.0f/g, dstSize);

  cudaSafeCall(cudaMemcpy(dst, d_dst, dst_mem, cudaMemcpyDeviceToHost));
  cudaSafeCall(cudaFree(d_src));
  cudaSafeCall(cudaFree(d_dst));

  return true;
}
} // namespace gpu

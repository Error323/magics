#include "init.h"
#include "../Verbose.hpp"

#include <stdio.h>
#include <stdlib.h>

namespace gpu
{

bool initCuda(void)
{
  int device;
  bool has_cuda = (cudaGetDevice(&device)) == cudaSuccess;
  if (!has_cuda)
    return false;

  cudaDeviceProp prop;
  cudaSafeCall(cudaGetDeviceProperties(&prop, device));

  DebugLine("");
  DebugLine("Initialized " << prop.name);
  DebugLine("  Compute capability:   " << prop.major << "." << prop.minor);
  DebugLine("  Multiprocessors:      " << prop.multiProcessorCount);
  DebugLine("  Global memory size:   " << prop.totalGlobalMem);
  DebugLine("  Constant memory size: " << prop.totalConstMem);
  DebugLine("  Concurrent kernels:   " << prop.concurrentKernels);
  DebugLine("  Shared mem per block: " << prop.sharedMemPerBlock);
  DebugLine("  Warpsize:             " << prop.warpSize);
  DebugLine("  Free/Total:           " << freeMemory() << "/" << totalMemory());
  DebugLine("");

  return true;
}

void error(const char *error_string, const char *file, const int line, const char *func)
{
  ErrorLine("\n\nCudaError: " << error_string);
  ErrorLine("  File:  " << file);
  ErrorLine("  Line:  " << line);
  ErrorLine("  Func:  " << func);
  exit(1);
}

size_t freeMemory()
{
  size_t free_memory, total_memory;
  cudaSafeCall(cudaMemGetInfo(&free_memory, &total_memory));
  return free_memory;
}

size_t totalMemory()
{
  size_t free_memory, total_memory;
  cudaSafeCall(cudaMemGetInfo(&free_memory, &total_memory));
  return total_memory;
}

}

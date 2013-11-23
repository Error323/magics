#ifndef INIT_CUH
#define INIT_CUH

#include <cuda_runtime_api.h>

#define cudaSafeCall(expr) gpu::___cudaSafeCall(expr, __FILE__, __LINE__, __func__)

namespace gpu
{
void error(const char *error_string, const char *file, const int line, const char *func = "");

static inline void ___cudaSafeCall(cudaError_t err, const char *file, const int line, const char *func = "")
{
    if (cudaSuccess != err)
        gpu::error(cudaGetErrorString(err), file, line, func);
}

bool initCuda(void);
size_t freeMemory(void);
size_t totalMemory(void);
}

#endif // INIT_CUH

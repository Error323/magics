#ifndef PROCESS_CUH
#define PROCESS_CUH

#include <cuda_runtime_api.h>
#include <curand.h>

#define NUM_PARENTS 20
#define BLOCK_DIM_1D 32
#define THREAD_DIM_1D 64
#define POOL_SIZE 2048
#define C64(x) x##ull

typedef unsigned long long U64;
typedef unsigned int U32;
typedef unsigned char U8;

namespace gpu
{
__device__ int Transform(U64 board, const U64 magic);
__global__ void InitPool(U64 *magics, U64 *randoms, int target_bits);
__global__ void ComputeFitness(U64 *magics, int *fitness, U64 *used_list, U64 *solution, int *sum, int n, int m);
__global__ void ComputeFitness(U64 *magics, int *fitness, U64 *used_list, U64 *solution, int *sum, int n, int m);
__global__ void SelectParents(U64 *magics, int *fitness, int *sum, U64 *parents, U32 *randoms);
__global__ void CreateOffspring(U64 *magics, U64 *parents, U64 *rand64);
bool process(int target_bits, int max_bits, const U64 *block, const U64 *attack);
}

#endif // PROCESS_CUH

#ifndef PROCESS_CUH
#define PROCESS_CUH

#include <cuda_runtime_api.h>
#include <curand.h>

#define NUM_ISLANDS 64
#define NUM_INDIVIDUALS 512
#define C64(x) x##ull

typedef unsigned long long U64;
typedef unsigned int U32;
typedef unsigned char U8;

namespace gpu
{
__device__ int Transform(U64 board, const U64 magic);
__global__ void InitPool(U64 *magics, U64 *randoms, int target_bits);
__global__ void SelectParents(U64 *magics, U64 *parents, U32 *collisions, U64 *used, int n, int m);
__global__ void CreateOffspring(U64 *magics, U64 *parents, U64 *rand64, int target_bits);
bool Process(int target_bits, int max_bits, const U64 *block, const U64 *attack, U64 &solution);
}

#endif // PROCESS_CUH

#include "process.h"
#include "init.h"
#include "../Verbose.hpp"
#include "../timer.h"

#include <math.h>
#include <math_constants.h>

// Algorithmic description:
// 1. initialize pool
//
// loop until fitness is good enough:
//   1. compute fitness for each chromosome
//   2. select N best parents
//   3. create new pool with children from N parents

extern bool stopped;
namespace gpu
{
__constant__ U64 block_list[1 << 12];
__constant__ U64 attack_list[1 << 12];

__device__ int Transform(U64 board, const U64 magic)
{
  board *= magic;
  board >>= (magic >> 58);
  return (int) board;
}

__global__ void InitPool(U64 *magics, U64 *randoms, int target_bits)
{
  int stride = blockIdx.x * blockDim.x;
  int id = stride + threadIdx.x;

  U64 magic = randoms[id] & randoms[id+stride] & randoms[id+stride*2] & C64(0x3ffffffffffffff);
  U64 shift = 64 - target_bits;
  magic |= shift << 58;
  magics[id] = magic;
}

__global__ void SelectParents(U64 *magics, U64 *parents, U32 *collisions, U64 *used, int n, int m)
{
  int id = blockIdx.x * blockDim.x + threadIdx.x;
  __shared__ U64 smagics[NUM_INDIVIDUALS];
  __shared__ U32 scollisions[NUM_INDIVIDUALS];
  U64 magic = magics[id];
  U32 col = 0;
  int start = id*m;
  int index;

  #pragma unroll 32
  for (int i = 0; i < n; i++)
  {
    index = Transform(block_list[i], magic);

    if (used[start + index] == C64(0))
      used[start + index] = attack_list[i];
    else
    if (used[start + index] != attack_list[i])
      col++;
  }

  scollisions[threadIdx.x] = col;
  smagics[threadIdx.x] = magic;
  __syncthreads();

  for (int offset = blockDim.x / 2; offset > 0; offset >>= 1)
  {
    if (threadIdx.x < offset)
    {
      if (scollisions[threadIdx.x] > scollisions[threadIdx.x + offset])
      {
        scollisions[threadIdx.x] = scollisions[threadIdx.x + offset];
        smagics[threadIdx.x] = smagics[threadIdx.x + offset];
      }
    }

    __syncthreads();
  }

  // thread 0 writes the final result
  if (threadIdx.x == 0)
  {
    collisions[blockIdx.x] = scollisions[0];
    parents[blockIdx.x] = smagics[0];
  }
}

__global__ void CreateOffspring(U64 *magics, U64 *parents, U64 *rand64, int target_bits)
{
  __shared__ U64 sparents[NUM_ISLANDS];

  if (threadIdx.x < NUM_ISLANDS)
    sparents[threadIdx.x] = parents[threadIdx.x];

  __syncthreads();

  int id = blockIdx.x * blockDim.x + threadIdx.x;
  U64 child;
  U64 r1     = rand64[id];
  U64 r2     = rand64[id+blockIdx.x*blockDim.x];
  U64 r3     = rand64[id+2*blockIdx.x*blockDim.x];
  int a = r1 % NUM_ISLANDS;
  int b = r2 % NUM_ISLANDS;
  int c = r3 % NUM_ISLANDS;

  // add random child
  if (a == c || b == c)
  {
    child = r1 & r2 & r3 & C64(0x3ffffffffffffff);
    U64 shift = 64 - target_bits;
    child |= shift << 58;
  }
  // have some good old binary sex
  else
  {
    U64 father = sparents[a];
    U64 mother = sparents[b];
    /*
    int crossover = (r3 % 62) + 1; 
    U64 father_side = (C64(1) << crossover) - 1;
    */
    U64 father_side = r1 ^ r2;
    child = (father & father_side) | (mother & ~father_side);
    child ^= (r1 & r2 & r3 & C64(0x3ffffffffffffff));
  }
  magics[id] = child;
}

bool Process(int target_bits, int max_bits, const U64 *block, const U64 *attack, U64 &solution)
{
  int n = 1 << max_bits;
  int m = 1 << target_bits;
  bool found = false;
  cudaSafeCall(cudaDeviceReset());

  // Store block and attack list in constant memory on device
  cudaSafeCall(cudaMemcpyToSymbol(block_list, block, sizeof(U64)*n));
  cudaSafeCall(cudaMemcpyToSymbol(attack_list, attack, sizeof(U64)*n));

  // Create random number generator
  curandGenerator_t rnd_gen;
  curandCreateGenerator(&rnd_gen, CURAND_RNG_PSEUDO_MTGP32);
  curandSetPseudoRandomGeneratorSeed(rnd_gen, time(0));

  // Allocate all required memory
  U64 *d_rand64;
  U64 *d_magics;
  U64 *d_used;
  U64 *d_parents;
  U32 *d_collisions;
  U64 h_parents[NUM_ISLANDS];
  U32 h_collisions[NUM_ISLANDS];
  cudaSafeCall(cudaMalloc((void**)&d_rand64, NUM_ISLANDS*NUM_INDIVIDUALS*3*sizeof(U64)));
  cudaSafeCall(cudaMalloc((void**)&d_magics, NUM_ISLANDS*NUM_INDIVIDUALS*sizeof(U64)));
  cudaSafeCall(cudaMalloc((void**)&d_collisions, NUM_ISLANDS*sizeof(U32)));
  cudaSafeCall(cudaMalloc((void**)&d_parents, NUM_ISLANDS*sizeof(U64)));
  cudaSafeCall(cudaMalloc((void**)&d_used, NUM_ISLANDS*NUM_INDIVIDUALS*m*sizeof(U64)));

  // Generate random numbers for the pool
  curandGenerate(rnd_gen, (U32*)d_rand64, NUM_ISLANDS*NUM_INDIVIDUALS*3*2);
  cudaSafeCall(cudaDeviceSynchronize());

  U32 generation = 0;
  U32 counter = 0;

  double start_time = timer::GetRealTime();
  char unit[4] = {'K','M','G','T'};

  // Initialize the pool
  InitPool<<<NUM_ISLANDS, NUM_INDIVIDUALS>>>(d_magics, d_rand64, target_bits);
  U32 collisions = 10000;
  U64 best = C64(0);
  while (!stopped && !found)
  {
    // Regenerate randoms
    curandGenerate(rnd_gen, (U32*)d_rand64, NUM_ISLANDS*NUM_INDIVIDUALS*3*2);
    cudaSafeCall(cudaDeviceSynchronize());

    double time = timer::GetRealTime() - start_time;
    if (time > 5)
    {
      double mps = counter / time;
      int u = floor(log10(mps)) / 3;
      u = std::max(std::min(u, 4), 1);
      mps /= pow(10, u*3);
      printf("G %d\tS %0.2f%c m/s C %u M 0x%llxull\n", generation, mps, unit[u-1], collisions, best);
      start_time += time;
      counter = 0;
    }

    // Compute fitness
    cudaSafeCall(cudaMemset(d_used, C64(0), NUM_ISLANDS*NUM_INDIVIDUALS*m*sizeof(U64)));
    cudaSafeCall(cudaMemset(d_collisions, 10000, NUM_ISLANDS*sizeof(U32)));
    SelectParents<<<NUM_ISLANDS, NUM_INDIVIDUALS>>>(d_magics, d_parents, d_collisions, d_used, n, m);
    
    // Check for solution
    cudaSafeCall(cudaMemcpy(h_parents, d_parents, NUM_ISLANDS*sizeof(U64), cudaMemcpyDeviceToHost));
    cudaSafeCall(cudaMemcpy(h_collisions, d_collisions, NUM_ISLANDS*sizeof(U32), cudaMemcpyDeviceToHost));

    for (int i = 0; i < NUM_ISLANDS; i++)
    {
      found = h_collisions[i] == 0;
      if (found)
      {
        solution = h_parents[i];
        break;
      }
      if (h_collisions[i] <= collisions)
      {
        collisions = h_collisions[i];
        best = h_parents[i];
      }
    }

    // Create offspring
    if (!found)
    {
      CreateOffspring<<<NUM_ISLANDS, NUM_INDIVIDUALS>>>(d_magics, d_parents, d_rand64, target_bits);
      generation++;
      counter += NUM_ISLANDS*NUM_INDIVIDUALS;
    }
  }
  
  // Free allocated cuda memory
  curandDestroyGenerator(rnd_gen);
  cudaSafeCall(cudaFree(d_rand64));
  cudaSafeCall(cudaFree(d_magics));
  cudaSafeCall(cudaFree(d_used));
  cudaSafeCall(cudaFree(d_collisions));
  cudaSafeCall(cudaFree(d_parents));

  return found;
}
} // namespace gpu

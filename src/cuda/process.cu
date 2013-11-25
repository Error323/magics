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
__constant__ U64 block_list[4096];
__constant__ U64 attack_list[4096];

__device__ int Transform(U64 board, const U64 magic)
{
  board *= magic;
  board >>= (magic >> 58);
  return (int) board;
}

__global__ void InitPool(U64 *magics, U64 *randoms, int target_bits)
{
  int id = blockIdx.x * THREAD_DIM_1D + threadIdx.x;

  if (id >= POOL_SIZE)
    return;

  U64 magic = randoms[id] & randoms[id+POOL_SIZE] & randoms[id+POOL_SIZE*2] & C64(0x3ffffffffffffff);
  U64 shift = 64 - target_bits;
  magic |= shift << 58;
  magics[id] = magic;
}

__global__ void ComputeFitness(U64 *magics, int *fitness, U64 *used_list, U64 *solution, int *sum, int n, int m)
{
  int id = blockIdx.x * THREAD_DIM_1D + threadIdx.x;

  if (id >= POOL_SIZE)
    return;

  U64 magic = magics[id];
  U64 used;
  int collisions = 0;
  int start = id*m;
  int index;
  int f;

  for (int i = 0; i < n; i++)
  {
    index = Transform(block_list[i], magic);
    used = used_list[start + index];

    if (used == C64(0))
      used_list[start + index] = attack_list[i];
    else
    if (used != attack_list[i])
      collisions++;
  }

  // we found a perfect solution
  if (collisions == 0)
    atomicExch(solution, magic);

  // compute the sum of all the fitness
  f = n - collisions;
  atomicAdd(sum, f);
  fitness[id] = f;
}

__global__ void SelectParents(U64 *magics, int *fitness, int *sum, U64 *parents, U32 *randoms)
{
  int pid = threadIdx.x;
  int r   = randoms[pid] % (*sum);
  int s   = fitness[0];
  for (int i = 1; i < POOL_SIZE; i++)
  {
    if (r < s)
    {
      parents[pid] = magics[i-1];
      break;
    }
    s += fitness[i];
  }
}

__global__ void CreateOffspring(U64 *magics, U64 *parents, U64 *rand64)
{
  int id = blockIdx.x * THREAD_DIM_1D + threadIdx.x;

  if (id >= POOL_SIZE)
    return;

  U64 child;
  U64 r1 = rand64[id];
  U64 r2 = rand64[id+POOL_SIZE];
  U64 r3 = rand64[id+2*POOL_SIZE];
  U64 father = parents[r1 % NUM_PARENTS];
  U64 mother = parents[r2 % NUM_PARENTS];
  int crossover = r3 % 64;
  U64 father_side = (C64(1) << crossover) - 1;
  child = (father & father_side) | (mother & ~father_side);
  child ^= (r1 & r2 & r3 & C64(0x3ffffffffffffff));
  magics[id] = child;
}

bool process(int target_bits, int max_bits, const U64 *block, const U64 *attack)
{
  int n = 1 << max_bits;
  int m = 1 << target_bits;

  // Store block and attack list in constant memory on device
  cudaSafeCall(cudaMemcpyToSymbol(block_list, block, sizeof(U64)*n));
  cudaSafeCall(cudaMemcpyToSymbol(attack_list, attack, sizeof(U64)*n));

  // Create random number generator
  curandGenerator_t rnd_gen;
  curandCreateGenerator(&rnd_gen, CURAND_RNG_PSEUDO_MTGP32);
  curandSetPseudoRandomGeneratorSeed(rnd_gen, time(0));

  // Generate random numbers for the pool
  U64 *d_rand64;
  cudaSafeCall(cudaMalloc((void**)&d_rand64, POOL_SIZE*3*sizeof(U64)));
  curandGenerate(rnd_gen, (U32*)d_rand64, POOL_SIZE*3*2);
  cudaSafeCall(cudaDeviceSynchronize());

  // Allocate all required memory
  U64 *d_magics;
  U64 *d_used, *d_solution;
  int *d_fitness, *d_sum;
  U64 *d_parents;
  U32 *d_rand32 = (U32*) d_rand64;
  cudaSafeCall(cudaMalloc((void**)&d_magics, POOL_SIZE*sizeof(U64)));
  cudaSafeCall(cudaMalloc((void**)&d_fitness, POOL_SIZE*sizeof(int)));
  cudaSafeCall(cudaMalloc((void**)&d_sum, sizeof(int)));
  cudaSafeCall(cudaMalloc((void**)&d_used, POOL_SIZE*m*sizeof(U64)));
  cudaSafeCall(cudaMalloc((void**)&d_solution, sizeof(U64)));
  cudaSafeCall(cudaMalloc((void**)&d_parents, NUM_PARENTS*sizeof(U64)));

  U64 solution = C64(0);
  U32 generation = 0;
  U32 counter = 0;

  double start_time = timer::GetRealTime();
  char unit[4] = {'K','M','G','T'};
  cudaSafeCall(cudaMemset(d_solution, C64(0), sizeof(U64)));

  // Initialize the pool
  InitPool<<<BLOCK_DIM_1D, THREAD_DIM_1D>>>(d_magics, d_rand64, target_bits);
  while (!stopped)
  {
    // Regenerate randoms
    curandGenerate(rnd_gen, (U32*)d_rand64, POOL_SIZE*3*2);

    double time = timer::GetRealTime() - start_time;
    if (time > 5)
    {
      double mps = counter / time;
      int u = floor(log10(mps)) / 3;
      u = std::max(std::min(u, 4), 1);
      mps /= pow(10, u*3);
      printf("G %d\tS %0.2f%c m/s\n", generation, mps, unit[u-1]);
      start_time += time;
      counter = 0;
    }

    // Compute fitness
    cudaSafeCall(cudaMemset(d_used, C64(0), POOL_SIZE*m*sizeof(U64)));
    cudaSafeCall(cudaMemset(d_sum, 0, sizeof(int)));
    ComputeFitness<<<BLOCK_DIM_1D, THREAD_DIM_1D>>>(d_magics, d_fitness, d_used, d_solution, d_sum, n, m);
    
    // Check for solution
    cudaSafeCall(cudaMemcpy(&solution, d_solution, sizeof(U64), cudaMemcpyDeviceToHost));
    if (solution != C64(0))
    {
      fprintf(stderr, "Solution found: 0x%llxull\n", solution);
      break;
    }


    // Select best <N> parents
    cudaSafeCall(cudaDeviceSynchronize());
    SelectParents<<<1, NUM_PARENTS>>>(d_magics, d_fitness, d_sum, d_parents, d_rand32);

    // Create offspring
    CreateOffspring<<<BLOCK_DIM_1D, THREAD_DIM_1D>>>(d_magics, d_parents, d_rand64);

    generation++;
    counter += POOL_SIZE;
  }

  
  // Free allocated cuda memory
  curandDestroyGenerator(rnd_gen);
  cudaSafeCall(cudaFree(d_rand64));
  cudaSafeCall(cudaFree(d_magics));
  cudaSafeCall(cudaFree(d_used));
  cudaSafeCall(cudaFree(d_fitness));
  cudaSafeCall(cudaFree(d_sum));
  cudaSafeCall(cudaFree(d_solution));
  cudaSafeCall(cudaFree(d_parents));
  cudaSafeCall(cudaDeviceReset());

  return true;
}
} // namespace gpu

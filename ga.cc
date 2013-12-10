#include <vector>
#include <random>

#include <immintrin.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <signal.h>
#include <stdio.h>
#include <getopt.h>
#include <assert.h>
#include <string.h>

#include "timer.h"
#include "align.h"

#define NUM_PARENTS 8
#define NUM_RANDOM 4
#define POPULATION_SIZE 1000

typedef unsigned long long U64;
typedef unsigned char U8;

#define C64(x) x##ull
#define RAND_FLT()    (dr(mt))
#define RAND_PARENT() (dp(mt))

// globals
std::mt19937_64 mt;
std::uniform_int_distribution<U64> d64;
std::uniform_int_distribution<U64> dp(0, NUM_PARENTS-1);
std::uniform_real_distribution<> dr(0, 1);
int square;
int is_bishop;
U64 mask;
U64 magic_seed;
unsigned int target_bits;
unsigned int max_bits;
unsigned int min_bits;
float fitness_sum;
std::vector<U64, AlignmentAllocator<U64, 32>> attack_list;
std::vector<U64, AlignmentAllocator<U64, 32>> block_list;
std::vector<U64, AlignmentAllocator<U64, 32>> used_list;

template<int N=1> U64 R64()
{
  U64 r = d64(mt);
  for (int i = 1; i < N; i++)
    r &= d64(mt);
  return r;
}

class Chromosome
{
public:
  Chromosome(): magic(C64(0)), collisions(0), fitness(0.0f) {}
  Chromosome(U64 m): magic(m), collisions(0), fitness(0.0f) {}
  Chromosome(U64 m, int f): magic(m), collisions(0), fitness(f) {}
  Chromosome(const Chromosome &c): magic(c.magic), collisions(c.collisions), fitness(c.fitness) {}

  Chromosome &operator=(const Chromosome &c)
  {
    if (this != &c)
    {
      magic = c.magic;
      fitness = c.fitness;
      collisions = c.collisions;
    }

    return *this;
  }

  U64 magic;
  int collisions;
  float fitness;
};

const int r_dist_attack_sets[64] =
{
  49, 42, 70,  84,  84,  70,  42, 49,
  42, 36, 60,  72,  72,  60,  36, 42,
  70, 60, 100, 120, 120, 100, 60, 70,
  84, 72, 120, 144, 144, 120, 72, 84,
  84, 72, 120, 144, 144, 120, 72, 84,
  70, 60, 100, 120, 120, 100, 60, 70,
  42, 36, 60,  72,  72,  60,  36, 42,
  49, 42, 70,  84,  84,  70,  42, 49
};

const int b_dist_attack_sets[64] =
{
  7,  6,  10, 12,  12,  10, 6,  7,
  6,  6,  10, 12,  12,  10, 6,  6,
  10, 10, 40, 48,  48,  40, 10, 10,
  12, 12, 48, 108, 108, 48, 12, 12,
  12, 12, 48, 108, 108, 48, 12, 12,
  10, 10, 40, 48,  48,  40, 10, 10,
  6,  6,  10, 12,  12,  10, 6,  6,
  7,  6,  10, 12,  12,  10, 6,  7
};

const int bit_table[64] =
{
  63, 30, 3 , 32, 25, 41, 22, 33,
  15, 50, 42, 13, 11, 53, 19, 34,
  61, 29, 2 , 51, 21, 43, 45, 10,
  18, 47, 1 , 54, 9 , 57, 0 , 35,
  62, 31, 40, 4 , 49, 5 , 52, 26,
  60, 6 , 23, 44, 46, 27, 56, 16,
  7 , 39, 48, 24, 59, 14, 12, 55,
  38, 28, 58, 20, 37, 17, 36, 8
};

void print(const U64 &inBoard)
{
  for (int i = 7; i >= 0; i--)
  {
    printf("%d", i+1);
    U8 line = (inBoard >> (i * 8)) & 0xff;

    for (int j = 0; j < 8; j++)
      printf(" %s", ((line >> j) & 1) == 1 ? "1" : ".");

    printf("\n");
  }

  printf("  a b c d e f g h\n");
}

U64 Magic(const U64 bits)
{
  U64 shift = 64 - bits;
  U64 magic = R64<3>() & C64(0x3ffffffffffffff);
  magic |= shift << 58;

  return magic;
}

int count_1s(U64 b)
{
  int r;
  for (r = 0; b; r++, b &= b - 1);
  return r;
}

int pop_1st_bit(U64 *bb)
{
  U64 b = *bb ^ (*bb - 1);
  unsigned int fold = (unsigned) ((b & 0xffffffff) ^ (b >> 32));
  *bb &= (*bb - 1);
  return bit_table[(fold * 0x783a9b23) >> 26];
}

U64 index_to_U64(int index, int bits, U64 m)
{
  int i, j;
  U64 result = C64(0);

  for (i = 0; i < bits; i++)
  {
    j = pop_1st_bit(&m);

    if (index & (1 << i))
      result |= (C64(1) << j);
  }

  return result;
}

U64 rmask(int sq)
{
  U64 result = C64(0);
  int rk = sq / 8, fl = sq % 8, r, f;

  for (r = rk + 1; r <= 6; r++)
    result |= (C64(1) << (fl + r * 8));

  for (r = rk - 1; r >= 1; r--)
    result |= (C64(1) << (fl + r * 8));

  for (f = fl + 1; f <= 6; f++)
    result |= (C64(1) << (f + rk * 8));

  for (f = fl - 1; f >= 1; f--)
    result |= (C64(1) << (f + rk * 8));

  return result;
}

U64 bmask(int sq)
{
  U64 result = C64(0);
  int rk = sq / 8, fl = sq % 8, r, f;

  for (r = rk + 1, f = fl + 1; r <= 6 && f <= 6; r++, f++)
    result |= (C64(1) << (f + r * 8));

  for (r = rk + 1, f = fl - 1; r <= 6 && f >= 1; r++, f--)
    result |= (C64(1) << (f + r * 8));

  for (r = rk - 1, f = fl + 1; r >= 1 && f <= 6; r--, f++)
    result |= (C64(1) << (f + r * 8));

  for (r = rk - 1, f = fl - 1; r >= 1 && f >= 1; r--, f--)
    result |= (C64(1) << (f + r * 8));

  return result;
}

U64 ratt(int sq, U64 block)
{
  U64 result = C64(0);
  int rk = sq / 8, fl = sq % 8, r, f;

  for (r = rk + 1; r <= 7; r++)
  {
    result |= (C64(1) << (fl + r * 8));

    if (block & (C64(1) << (fl + r * 8)))
      break;
  }

  for (r = rk - 1; r >= 0; r--)
  {
    result |= (C64(1) << (fl + r * 8));

    if (block & (C64(1) << (fl + r * 8)))
      break;
  }

  for (f = fl + 1; f <= 7; f++)
  {
    result |= (C64(1) << (f + rk * 8));

    if (block & (C64(1) << (f + rk * 8)))
      break;
  }

  for (f = fl - 1; f >= 0; f--)
  {
    result |= (C64(1) << (f + rk * 8));

    if (block & (C64(1) << (f + rk * 8)))
      break;
  }

  return result;
}

U64 batt(int sq, U64 block)
{
  U64 result = C64(0);
  int rk = sq / 8, fl = sq % 8, r, f;

  for (r = rk + 1, f = fl + 1; r <= 7 && f <= 7; r++, f++)
  {
    result |= (C64(1) << (f + r * 8));

    if (block & (C64(1) << (f + r * 8)))
      break;
  }

  for (r = rk + 1, f = fl - 1; r <= 7 && f >= 0; r++, f--)
  {
    result |= (C64(1) << (f + r * 8));

    if (block & (C64(1) << (f + r * 8)))
      break;
  }

  for (r = rk - 1, f = fl + 1; r >= 0 && f <= 7; r--, f++)
  {
    result |= (C64(1) << (f + r * 8));

    if (block & (C64(1) << (f + r * 8)))
      break;
  }

  for (r = rk - 1, f = fl - 1; r >= 0 && f >= 0; r--, f--)
  {
    result |= (C64(1) << (f + r * 8));

    if (block & (C64(1) << (f + r * 8)))
      break;
  }

  return result;
}


int transform(U64 board, const U64 magic)
{
  board *= magic;
  board >>= (magic >> 58);

  return static_cast<int>(board);
}

void ComputeFitness(Chromosome &chromosome)
{
  #if defined __AVX__
  __m256i *dst = reinterpret_cast<__m256i*>(used_list.data());
  for (int i = 0, n = used_list.size()/4; i < n; i++)
    dst[i] = _mm256_set1_epi64x(0);
  #elif defined __SSE4_1__
  __m128i *dst = reinterpret_cast<__m128i*>(used_list.data());
  for (int i = 0, n = used_list.size()/2; i < n; i++)
    dst[i] = _mm_set1_epi64x(0);
  #else
  used_list.assign(used_list.size(), C64(0));
  #endif 
  int index, n, i;
  chromosome.collisions = 0;

  n = (1 << max_bits);
  for (i = 0; i < n; i++)
  {
    index = transform(block_list[i], chromosome.magic);

    if (used_list[index] == C64(0))
      used_list[index] = attack_list[i];
    else
    if (used_list[index] != attack_list[i])
      chromosome.collisions++;
  }

  chromosome.fitness = n - chromosome.collisions;
}

void InitializePopulation(std::vector<Chromosome> &pool)
{
  for (int i = 0; i < POPULATION_SIZE; i++)
  {
    pool[i].magic = Magic(target_bits);
    ComputeFitness(pool[i]);
  }

  if (magic_seed != C64(0))
  {
    U64 shift = 64-target_bits;
    if (magic_seed >> 58 != shift)
    {
      magic_seed &= C64(0x3ffffffffffffff);
      magic_seed |= shift << 58;
    }
    
    pool[0].magic = magic_seed;
    ComputeFitness(pool[0]);
  }
}

void GetBestSolution(std::vector<Chromosome> &pool, Chromosome &solution)
{
  fitness_sum = 0.0f;
  float best = 0;
  for (int i = 0; i < POPULATION_SIZE; i++)
  {
    fitness_sum += pool[i].fitness;
    if (pool[i].fitness > best)
    {
      best = pool[i].fitness;
      solution = pool[i];
    }
  }
}

void SelectParents(const std::vector<Chromosome> &pool, std::vector<Chromosome> &parents, const Chromosome &solution)
{
  // Add parents using roulette wheel selection
  int num_parents = NUM_PARENTS-NUM_RANDOM;
  int randoms[num_parents];
  bool done[num_parents];
  for (int i = 1; i < num_parents; i++)
  {
    randoms[i] = RAND_FLT() * fitness_sum;
    done[i] = false;
  }

  // Hard add the previous best solution
  parents[0] = solution;
  done[0] = true;

  int sum = 0;
  bool can_stop = true;
  for (int i = 0; i < POPULATION_SIZE; i++)
  {
    sum += pool[i].fitness;
    for (int j = 1; j < num_parents; j++)
    {
      can_stop = can_stop && done[j];

      if (done[j])
        continue;

      if (randoms[j] < sum)
      {
        parents[j] = pool[i];
        done[j] = true;
      }
    }

    if (can_stop)
      break;
  }

  // Add some random parents
  for (int i = num_parents; i < NUM_PARENTS; i++)
  {
    parents[i].magic = Magic(target_bits);
    ComputeFitness(parents[i]);
  }
}

void GenerateOffspring(std::vector<Chromosome> &pool, const std::vector<Chromosome> &parents)
{
  // Put first <n> best parents in new pool
  for (int i = 0; i < NUM_PARENTS - NUM_RANDOM; i++)
    pool[i] = parents[i];

  // Create offspring
  for (int i = NUM_PARENTS - NUM_RANDOM; i < POPULATION_SIZE; i++)
  {
    // Select parents
    Chromosome &child = pool[i];
    int r1 = RAND_PARENT();
    int r2 = RAND_PARENT();
    while (r1 == r2)
      r2 = RAND_PARENT();
    const Chromosome &father = parents[r1];
    const Chromosome &mother = parents[r2];

    // Mate between father and mother
    U64 father_side = R64();
    child = (father.magic & father_side) | (mother.magic & ~father_side);

    // Mutate some bits
    child.magic ^= R64<3>() & C64(0x3ffffffffffffff);
    
    // Compute new fitness
    ComputeFitness(child);
  }
}

bool stopped;
void stop(int)
{
  stopped = true;
}

void print_and_exit(int ret)
{
  printf("Usage: magics [OPTION]...\n");
  printf(" Find magics for a given square.\n");
  printf(" Recover shift with (magic >> 58).\n");
  printf(" Prints solution to stderr.\n");
  printf(" Example: magics -s 0 -t 11\n\n");
  printf(" -h\tdisplay this help message\n");
  printf(" -s\tsquare to look for in {0,...,63}\n");
  printf(" -1\tlook for magics using max_bits-1\n");
  if (min_bits == 0 || max_bits == 0)
    printf(" -t\tbits to use\n");
  else
    printf(" -t\tbits to use in {%d,...,%d}\n", min_bits, max_bits);
  printf(" -b\tsearch for bishop, otherwise for rook\n");
  printf(" -g\tmax generations\n");
  printf(" -m\tseed population with this magic number\n");

  exit(ret);
}

int main(int argc, char **argv)
{
  std::random_device rd;
  mt.seed(rd());
  signal(SIGINT, stop);
  signal(SIGTERM, stop);

  target_bits = 0;
  square = 0;
  is_bishop = 0;
  min_bits = 0;
  max_bits = 0;
  magic_seed = C64(0);
  bool max_bits_min_one = false;
  int max_generations = 0;

  int c;
  while ((c = getopt(argc, argv, "1s:t:g:bhm:")) != -1)
  {
    switch (c)
    {
    case 's': square = atoi(optarg); break;
    case 't': target_bits = atoi(optarg); break;
    case 'm': magic_seed = strtoull(optarg, 0, 0); break;
    case 'b': is_bishop = 1; break;
    case '1': max_bits_min_one = true; break;
    case 'g': max_generations = atoi(optarg); break;
    case 'h': print_and_exit(EXIT_SUCCESS);
    case '?':
    default: print_and_exit(EXIT_FAILURE);
    }
  }

  if (argc == 1)
    print_and_exit(EXIT_FAILURE);

  if (square < 0 || square > 63)
  {
    fprintf(stderr, "Error: target square %d not in {0,...,63}\n", square);
    return EXIT_FAILURE;
  }

  mask = is_bishop ? bmask(square) : rmask(square);
  max_bits = count_1s(mask);
  min_bits = is_bishop ? b_dist_attack_sets[square] : r_dist_attack_sets[square];
  min_bits = ceil(log2(min_bits));

  if (target_bits == 0) target_bits = max_bits;
  if (max_bits_min_one) target_bits = max_bits-1;

  if (target_bits < min_bits || target_bits > max_bits)
  {
    fprintf(stderr, "Error: target bits %d not in {%d,...,%d}\n", target_bits, min_bits, max_bits);
    return EXIT_FAILURE;
  }

  attack_list.resize(1 << max_bits);
  block_list.resize(1 << max_bits);
  used_list.resize(1 << target_bits);

  for (int i = 0; i < (1 << max_bits); i++)
  {
    block_list[i] = index_to_U64(i, max_bits, mask);
    attack_list[i] = is_bishop ? batt(square, block_list[i]) : 
                     ratt(square, block_list[i]);
  }

  std::vector<Chromosome> pool(POPULATION_SIZE);
  std::vector<Chromosome> parents(NUM_PARENTS);
  Chromosome solution;
  InitializePopulation(pool);

  printf("Generating magic for '%s' on square %c%d using %d <= (%d) <= %d bits\n", 
          is_bishop ? "bishop" : "rook", char(square%8+65),
          square/8+1, min_bits, target_bits, max_bits);

  char cmd_line[256];
  for (int i = 0, j = 0; i < argc; i++)
  {
    strcpy(&cmd_line[j], argv[i]);
    j += strlen(argv[i]);
  }

  constexpr char unit[] = {'K', 'M', 'G', 'T'};
  int generation = 0;
  stopped = false;
  bool solution_found = false;
  double start_time = GetRealTime();
  int counter = 0;
  while (!stopped)
  {
    GetBestSolution(pool, solution);

    //ComputeCollisions(solution);
    solution_found = solution.collisions == 0;

    counter++;
    double time = GetRealTime() - start_time;
    if (time >= 5.0 || solution_found)
    {
      double mps = (POPULATION_SIZE * counter) / time;
      int u = floor(log10(mps)) / 3;
      u = std::max(std::min(u, 4), 1);
      mps /= pow(10, u*3);
      printf("G %10d\tS %0.2f%c C %d\n", generation, mps, unit[u-1], solution.collisions);
      start_time += time;
      counter = 0;
      fflush(stdout);
    }

    if (solution_found)
      break;

    SelectParents(pool, parents, solution);
    GenerateOffspring(pool, parents);
    generation++;

    if (generation > max_generations && max_generations > 0)
      break;
  }

  if (solution_found)
  {
    assert(target_bits == 64-(solution.magic >> 58));
    fprintf(stderr, "G %d\t0x%llxull\t%s\t %c%d\t%d\n",
            generation, solution.magic, is_bishop ? "bishop" : "rook",
            char(square%8+65), square/8+1, square);
  }
  else
    fprintf(stderr, "FAIL\t%s\t %c%d\t%d\n",
            is_bishop ? "bishop" : "rook",
            char(square%8+65), square/8+1, square);

  return EXIT_SUCCESS;
}

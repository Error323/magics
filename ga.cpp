#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include <algorithm>

#include <ctime>
#include <cmath>
#include <cstdlib>
#include <csignal>
#include <cstdio>

#define NUM_PARENTS 10
#define NUM_RANDOM 2
#define CROSSOVER_RATE 0.99f
#define MUTATION_RATE 0.01f
#define POPULATION_SIZE 1000
#define CHROMO_LENGTH 64

typedef unsigned long long U64;
typedef unsigned char U8;

#define C64(x) x##ULL
#define RAND_FLT()    (rand()/static_cast<float>(RAND_MAX))
#define RAND_INT(a,b) ((rand() % ((b)-(a) + 1)) + (a))

// globals
int square;
int is_bishop;
U64 mask;
int target_bits;
int max_bits;
int min_bits;
std::vector<U64> attack_list;
std::vector<U64> block_list;
std::vector<U64> used_list;

const int bit_table[64] = {
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
  std::stringstream ss; 
  for (int i = 7; i >= 0; i--)
  {
    ss << i + 1;
    U8 line = (inBoard >> (i * 8)) & 0xff;
    for (int j = 0; j < 8; j++)
      ss << (((line >> j) & 1) == 1 ? " 1" : " .");
    ss << "\n";
  }
  ss << "  a b c d e f g h\n";
  std::cout << ss.str();
}

U64 R64() 
{
  U64 u1, u2, u3, u4;
  u1 = (U64)(random()) & 0xFFFF; u2 = (U64)(random()) & 0xFFFF;
  u3 = (U64)(random()) & 0xFFFF; u4 = (U64)(random()) & 0xFFFF;
  return u1 | (u2 << 16) | (u3 << 32) | (u4 << 48);
}

U64 R64Few() 
{
  return R64() & R64() & R64();
}

int count_1s(U64 b) {
  int r;
  for(r = 0; b; r++, b &= b - 1);
  return r;
}

int pop_1st_bit(U64 *bb) {
  U64 b = *bb ^ (*bb - 1);
  unsigned int fold = (unsigned) ((b & 0xffffffff) ^ (b >> 32));
  *bb &= (*bb - 1);
  return bit_table[(fold * 0x783a9b23) >> 26];
}

U64 index_to_U64(int index, int bits, U64 m) {
  int i, j;
  U64 result = C64(0);
  for(i = 0; i < bits; i++) {
    j = pop_1st_bit(&m);
    if(index & (1 << i)) result |= (C64(1) << j);
  }
  return result;
}

U64 rmask(int sq) {
  U64 result = C64(0);
  int rk = sq/8, fl = sq%8, r, f;
  for(r = rk+1; r <= 6; r++) result |= (C64(1) << (fl + r*8));
  for(r = rk-1; r >= 1; r--) result |= (C64(1) << (fl + r*8));
  for(f = fl+1; f <= 6; f++) result |= (C64(1) << (f + rk*8));
  for(f = fl-1; f >= 1; f--) result |= (C64(1) << (f + rk*8));
  return result;
}

U64 bmask(int sq) {
  U64 result = C64(0);
  int rk = sq/8, fl = sq%8, r, f;
  for(r=rk+1, f=fl+1; r<=6 && f<=6; r++, f++) result |= (C64(1) << (f + r*8));
  for(r=rk+1, f=fl-1; r<=6 && f>=1; r++, f--) result |= (C64(1) << (f + r*8));
  for(r=rk-1, f=fl+1; r>=1 && f<=6; r--, f++) result |= (C64(1) << (f + r*8));
  for(r=rk-1, f=fl-1; r>=1 && f>=1; r--, f--) result |= (C64(1) << (f + r*8));
  return result;
}

U64 ratt(int sq, U64 block) {
  U64 result = C64(0);
  int rk = sq/8, fl = sq%8, r, f;
  for(r = rk+1; r <= 7; r++) {
    result |= (C64(1) << (fl + r*8));
    if(block & (C64(1) << (fl + r*8))) break;
  }
  for(r = rk-1; r >= 0; r--) {
    result |= (C64(1) << (fl + r*8));
    if(block & (C64(1) << (fl + r*8))) break;
  }
  for(f = fl+1; f <= 7; f++) {
    result |= (C64(1) << (f + rk*8));
    if(block & (C64(1) << (f + rk*8))) break;
  }
  for(f = fl-1; f >= 0; f--) {
    result |= (C64(1) << (f + rk*8));
    if(block & (C64(1) << (f + rk*8))) break;
  }
  return result;
}

U64 batt(int sq, U64 block) {
  U64 result = C64(0);
  int rk = sq/8, fl = sq%8, r, f;
  for(r = rk+1, f = fl+1; r <= 7 && f <= 7; r++, f++) {
    result |= (C64(1) << (f + r*8));
    if(block & (C64(1) << (f + r * 8))) break;
  }
  for(r = rk+1, f = fl-1; r <= 7 && f >= 0; r++, f--) {
    result |= (C64(1) << (f + r*8));
    if(block & (C64(1) << (f + r * 8))) break;
  }
  for(r = rk-1, f = fl+1; r >= 0 && f <= 7; r--, f++) {
    result |= (C64(1) << (f + r*8));
    if(block & (C64(1) << (f + r * 8))) break;
  }
  for(r = rk-1, f = fl-1; r >= 0 && f >= 0; r--, f--) {
    result |= (C64(1) << (f + r*8));
    if(block & (C64(1) << (f + r * 8))) break;
  }
  return result;
}


int transform(U64 b, U64 magic, int bits) 
{
  return (int)((b * magic) >> (64 - bits));
}

void InitializePopulation(std::vector<U64> &pool)
{
  for (int i = 0; i < POPULATION_SIZE; i++)
  {
    pool[i] = R64Few();
    if(count_1s((mask * pool[i]) & 0xFF00000000000000ULL) < min_bits)
      i--;
  }
}

int GetFitness(const U64 chromosome)
{
  used_list.assign(used_list.size(), C64(0));
  int bad_collisions = 0;
  int index;
  U64 attack;

  for (int i = 0; i < (1<<max_bits); i++)
  {
    attack = attack_list[i];
    index = transform(block_list[i], chromosome, target_bits);
    if (used_list[index] == C64(0))
      used_list[index] = attack;
    else
    if (used_list[index] != attack)
      bad_collisions++;
  }

  return (1<<max_bits) - bad_collisions;
}

bool SortChromosome(const U64 &a, const U64 &b)
{
  return GetFitness(a) > GetFitness(b);
}

void GetBestSolution(std::vector<U64> &pool, U64 &solution)
{
  std::sort(pool.begin(), pool.end(), SortChromosome);
  solution = pool[0];
}

void SelectParents(const std::vector<U64> &pool, std::vector<U64> &parents)
{
  //NOTE: Assuming sorted pool - descending order
  for (int i = 0; i < NUM_PARENTS-NUM_RANDOM; i++)
    parents[i] = pool[i];

  // Add some random parents
  for (int i = NUM_PARENTS-NUM_RANDOM; i < NUM_PARENTS; i++)
  {
    parents[i] = R64Few();
    if(count_1s((mask * parents[i]) & 0xFF00000000000000ULL) < min_bits)
      i--;
  }
}

void GenerateOffspring(std::vector<U64> &pool, const std::vector<U64> &parents)
{
  // Put best parents in new pool
  for (int i = 0; i < NUM_PARENTS-NUM_RANDOM; i++)
    pool[i] = parents[i];

  // Create offspring
  for (int i = NUM_PARENTS-NUM_RANDOM; i < POPULATION_SIZE; i++)
  {
    U64 &child = pool[i];
    const U64 &father = parents[RAND_INT(0, NUM_PARENTS-1)];
    const U64 &mother = parents[RAND_INT(0, NUM_PARENTS-1)];

    // Mate between father and mother
    if (RAND_FLT() < CROSSOVER_RATE)
    {
      int crossover = RAND_INT(0, CHROMO_LENGTH-1);
      U64 father_side = (C64(1) << crossover) - 1;
      child = (father&father_side) | (mother&~father_side);
    }

    // Create possible mutations
    for (int i = 0; i < CHROMO_LENGTH; i++)
      if (RAND_FLT() < MUTATION_RATE)
        child ^= (C64(1) << i);
  }
}

bool stopped;
void stop(int)
{
  stopped = true;
}

int main(int argc, char **argv)
{
  srand(time(0));
  signal(SIGINT, stop);

  min_bits = 7;
  target_bits = 11;
  square = 0;
  is_bishop = 0;

  switch (argc)
  {
    case 4:
      square = atoi(argv[1]);
      target_bits = atoi(argv[2]);
      is_bishop = atoi(argv[3]);
    break;
    case 3:
      square = atoi(argv[1]);
      target_bits = atoi(argv[2]);
    break;
    default:
      printf("Usage: magics <square> <target-bits> [is_bishop]\n");
      return EXIT_FAILURE;
    break;
  }

  mask = is_bishop ? bmask(square) : rmask(square);
  max_bits = count_1s(mask);

  attack_list.resize(1<<max_bits);
  block_list.resize(1<<max_bits);
  used_list.resize(1<<target_bits);

  for (int i = 0; i < (1<<max_bits); i++)
  {
    block_list[i] = index_to_U64(i, max_bits, mask);
    attack_list[i] = is_bishop ? batt(square, block_list[i]) : ratt(square, block_list[i]);
  }

  std::vector<U64> pool(POPULATION_SIZE);
  std::vector<U64> parents(NUM_PARENTS);
  U64 solution = C64(0);
  InitializePopulation(pool);

  int generation = 0, fitness;
  stopped = false;
  printf("Generating magic for '%s' using %d/%d bits\n", is_bishop?"bishop":"rook", target_bits, max_bits);
  print(C64(1) << square);
  while (!stopped)
  {
    GetBestSolution(pool, solution);
    fitness = GetFitness(solution);

    if (generation % 100 == 0 || fitness == (1<<max_bits))
      printf("Generation[%4d] bad collisions(%d): 0x%llxULL\n", generation, (1<<max_bits)-fitness, solution);

    if (fitness == (1<<max_bits))
    {
      printf("Found magic for '%s' using %d/%d bits\n", is_bishop?"bishop":"rook", target_bits, max_bits);
      print(C64(1) << square);
      break;
    }

    SelectParents(pool, parents);
    GenerateOffspring(pool, parents);
    generation++;
  }

  return EXIT_SUCCESS;
}

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
#define MUTATION_RATE 0.1f
#define POPULATION_SIZE 10000
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

class Chromosome
{
public:
  Chromosome(): magic(C64(0)), fitness(-1) {}
  Chromosome(U64 m): magic(m), fitness(-1) {}
  Chromosome(U64 m, int f): magic(m), fitness(f) {}
  Chromosome(const Chromosome &c): magic(c.magic), fitness(c.fitness) {}

  Chromosome &operator=(const Chromosome &c)
  {
    if (this != &c)
    {
      magic = c.magic;
      fitness = c.fitness;
    }

    return *this;
  }

  U64 magic;
  int fitness;
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
  u1 = (U64)(random()) & 0xFFFF;
  u2 = (U64)(random()) & 0xFFFF;
  u3 = (U64)(random()) & 0xFFFF;
  u4 = (U64)(random()) & 0xFFFF;
  return u1 | (u2 << 16) | (u3 << 32) | (u4 << 48);
}

U64 R64Few()
{
  return R64() & R64() & R64();
}

int count_1s(U64 b)
{
  int r;

  for (r = 0; b; r++, b &= b - 1)
    ;

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


int transform(U64 b, U64 magic, int bits)
{
  return (int)((b * magic) >> (64 - bits));
}

int GetFitness(const U64 chromosome)
{
  used_list.assign(used_list.size(), C64(0));
  int bad_collisions = 0;
  int index;
  U64 attack;

  for (int i = 0; i < (1 << max_bits); i++)
  {
    attack = attack_list[i];
    index = transform(block_list[i], chromosome, target_bits);

    if (used_list[index] == C64(0))
      used_list[index] = attack;
    else
      if (used_list[index] != attack)
        bad_collisions++;
  }

  return (1 << max_bits) - bad_collisions;
}

void InitializePopulation(std::vector<Chromosome> &pool)
{
  for (int i = 0; i < POPULATION_SIZE; i++)
  {
    pool[i].magic = R64Few();

    if (count_1s((mask * pool[i].magic) & 0xFF00000000000000ULL) < min_bits)
      i--;
    else
      pool[i].fitness = GetFitness(pool[i].magic);
  }
}

bool ChromosomeSorter(const Chromosome &a, const Chromosome &b)
{
  return a.fitness > b.fitness;
}

void GetBestSolution(std::vector<Chromosome> &pool, Chromosome &solution)
{
  std::sort(pool.begin(), pool.end(), ChromosomeSorter);
  solution = pool[0];
}

void SelectParents(const std::vector<Chromosome> &pool, std::vector<Chromosome> &parents)
{
  //NOTE: Assuming sorted pool - descending order
  for (int i = 0; i < NUM_PARENTS - NUM_RANDOM; i++)
    parents[i] = pool[i];

  // Add some random parents
  for (int i = NUM_PARENTS - NUM_RANDOM; i < NUM_PARENTS; i++)
  {
    parents[i].magic = R64Few();

    if (count_1s((mask * parents[i].magic) & 0xFF00000000000000ULL) < min_bits)
      i--;
    else
      parents[i].fitness = GetFitness(parents[i].magic);
  }
}

void GenerateOffspring(std::vector<Chromosome> &pool, const std::vector<Chromosome> &parents)
{
  // Put best parents in new pool
  for (int i = 0; i < NUM_PARENTS - NUM_RANDOM; i++)
    pool[i] = parents[i];

  // Create offspring
  for (int i = NUM_PARENTS - NUM_RANDOM; i < POPULATION_SIZE; i++)
  {
    Chromosome &child = pool[i];
    const Chromosome &father = parents[RAND_INT(0, NUM_PARENTS - 1)];
    const Chromosome &mother = parents[RAND_INT(0, NUM_PARENTS - 1)];

    // Mate between father and mother
    if (RAND_FLT() < CROSSOVER_RATE)
    {
      int crossover = RAND_INT(0, CHROMO_LENGTH - 1);
      U64 father_side = (C64(1) << crossover) - 1;
      child = (father.magic & father_side) | (mother.magic&~father_side);
    }

    // Create possible mutations
    for (int j = 0; j < CHROMO_LENGTH; j++)
      if (RAND_FLT() < MUTATION_RATE)
        child.magic ^= (C64(1) << j);

    child.fitness = GetFitness(child.magic);
  }
}

bool stopped;
void stop(int)
{
  stopped = true;
}

void print_and_exit()
{
  printf("Usage: magics <square> <target-bits> [is_bishop]\n");
  printf("       <square> in {0,...,63}\n");
  printf("       <target-bits> number of bits to look for magics\n");
  printf("       <is_bishop> in {1,0} bishop or rook\n");
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
  srand(time(0));
  signal(SIGINT, stop);
  signal(SIGTERM, stop);

  min_bits = 3;
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
    print_and_exit();
    break;
  }

  mask = is_bishop ? bmask(square) : rmask(square);
  max_bits = count_1s(mask);

  if (square < 0 || square > 63)
    print_and_exit();

  if (target_bits > max_bits)
    print_and_exit();

  attack_list.resize(1 << max_bits);
  block_list.resize(1 << max_bits);
  used_list.resize(1 << target_bits);

  for (int i = 0; i < (1 << max_bits); i++)
  {
    block_list[i] = index_to_U64(i, max_bits, mask);
    attack_list[i] = is_bishop ? batt(square, block_list[i]) : ratt(square, block_list[i]);
  }

  std::vector<Chromosome> pool(POPULATION_SIZE);
  std::vector<Chromosome> parents(NUM_PARENTS);
  Chromosome solution;
  InitializePopulation(pool);

  int generation = 0;
  stopped = false;
  printf("Generating magic for '%s' using %d/%d bits\n", is_bishop ? "bishop" : "rook", target_bits, max_bits);
  print(is_bishop ? batt(square, C64(0)) : ratt(square, C64(0)));
  bool solution_found = false;

  while (!stopped)
  {
    GetBestSolution(pool, solution);

    solution_found = solution.fitness == (1 << max_bits);

    if (solution_found)
      break;

    if (generation % 100 == 0)
      printf("Generation[%8d] bad collisions(%d): 0x%llxULL\n", generation, (1 << max_bits) - solution.fitness, solution.magic);

    SelectParents(pool, parents);
    GenerateOffspring(pool, parents);
    generation++;
  }

  if (solution_found)
    printf("Solution after %d generations: 0x%llxULL for '%s' square %d %d/%d bits.\n",
           generation, solution.magic, is_bishop ? "bishop" : "rook", square, target_bits, max_bits);
  else
    printf("No solution found after %d generations\n", generation);

  return EXIT_SUCCESS;
}

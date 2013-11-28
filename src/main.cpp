#include <vector>

#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <math.h>
#include <getopt.h>
#include <assert.h>
#include <string.h>
#include <x86intrin.h>

#include "timer.h"
#include "Verbose.hpp"
#include "cuda/init.h"
#include "cuda/process.h"

typedef unsigned long long U64;
typedef unsigned char U8;

#define C64(x) x##ull

// globals
int square;
int is_bishop;
U64 mask;
unsigned int target_bits;
unsigned int max_bits;
unsigned int min_bits;
std::vector<U64> attack_list;
std::vector<U64> block_list;

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

int transform(U64 board, const U64 magic)
{
  board *= magic;
  board >>= (magic >> 58);

  return static_cast<int>(board);
}

int collisions(U64 magic)
{
  int m = 1 << target_bits;
  U64 used_list[m];
  memset(used_list, 0, sizeof(U64)*m);
  int index, n, c = 0;

  n = (1 << max_bits);
  for (int i = 0; i < n; i++)
  {
    index = transform(block_list[i], magic);

    if (used_list[index] == C64(0))
      used_list[index] = attack_list[i];
    else
    if (used_list[index] != attack_list[i])
      c++;
  }

  return c;
}

int count_1s(register const U64 b)
{
  return __popcntq(b);
}

int pop_1st_bit(register U64 &b)
{
  int index = __bsfq(b);
  b &= b - 1;
  return index;
}

int idx_1st_bit(register const U64 b)
{
  return __bsfq(b);
}

U64 index_to_U64(int index, int bits, U64 m)
{
  int i, j;
  U64 result = C64(0);

  for (i = 0; i < bits; i++)
  {
    j = pop_1st_bit(m);

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

  exit(ret);
}

int main(int argc, char **argv)
{
  signal(SIGINT, stop);
  signal(SIGTERM, stop);

  target_bits = 0;
  square = 0;
  is_bishop = 0;
  min_bits = 0;
  max_bits = 0;
  bool max_bits_min_one = false;

  int c;
  while ((c = getopt(argc, argv, "1s:t:bh")) != -1)
  {
    switch (c)
    {
    case 's': square = atoi(optarg); break;
    case 't': target_bits = atoi(optarg); break;
    case 'b': is_bishop = 1; break;
    case '1': max_bits_min_one = true; break;
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

  Verbose::SetVerbosity(Verbose::DBG);
  gpu::initCuda();

  attack_list.resize(1 << max_bits);
  block_list.resize(1 << max_bits);

  for (int i = 0; i < (1 << max_bits); i++)
  {
    block_list[i] = index_to_U64(i, max_bits, mask);
    attack_list[i] = is_bishop ? batt(square, block_list[i]) : ratt(square, block_list[i]);
  }

  printf("Generating magic for '%s' on square %c%d using %d <= (%d) <= %d bits\n", 
          is_bishop ? "bishop" : "rook", char(square%8+65), square/8+1,
          min_bits, target_bits, max_bits);

  U64 solution;
  if (gpu::Process(target_bits, max_bits, &block_list[0], &attack_list[0], solution))
  {
    assert(target_bits == 64-(solution >> 58));
    assert(collisions(solution) == 0);
    fprintf(stderr, "0x%llxull\t%s\t %c%d\t%d\n",
            solution, is_bishop ? "bishop" : "rook",
            char(square%8+65), square/8+1, square);
  }
  else
  {
    fprintf(stderr, "FAIL\t%s\t %c%d\t%d\n",
            is_bishop ? "bishop" : "rook",
            char(square%8+65), square/8+1, square);
  }

  return EXIT_SUCCESS;
}

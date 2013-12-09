#ifndef BITBOARD_H
#define BITBOARD_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <x86intrin.h>

#include <vector>

using namespace std;

typedef unsigned long long U64;
typedef unsigned long U32;
typedef unsigned short U16;
typedef unsigned char U8;
typedef signed long long I64;
typedef signed long I32;
typedef signed short I16;
typedef signed char I8;

#define C64(x) x##ull

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

/// Pops the least significant bit of the board and returns index
inline int PopLSB(register U64 &b)
{
  int index = __bsfq(b);
  b &= b - 1;
  return index;
}

/// Find the index to the least significant bit
inline int LSBIndex(register U64 b)
{
  return __bsfq(b);
}

/// Count the number of set bits on a bitboard
inline int CountBits(register U64 b)
{
  return __popcntq(b);
}

U64 Index2U64(int index, int bits, U64 m);
U64 bmask(int sq);
U64 rmask(int sq);
U64 ratt(int sq, U64 block);
U64 batt(int sq, U64 block);

/// Print a bitboard to stdout
void Print1(const U64 board);

/// Print 2 bitboards to stdout beside eachother
void Print2(const U64 a, const U64 b);

/// Print 3 bitboards to stdout beside eachother
void Print3(const U64 a, const U64 b, const U64 c);

#endif // BITBOARD_H

#ifndef ISLAND_H
#define ISLAND_H

#include <vector>
#include <random>
#include <thread>

#include "bitboard.h"

using namespace std;

typedef pair<int, U64> chromosome;

class Island
{
public:
  static const int num_individuals = 1000;

  Island(int tbits, int nparents);
  ~Island();

  static void Init(const int maxbits, const bool bishop, const int square);
  chromosome RandChromosome();
  chromosome Populate();
  chromosome NextGen(const vector<chromosome> &parents, U64 best, int best_fitness);

private:
  int mTargetBits;
  vector<U64> mUsed;

  mt19937_64 mMt;
  uniform_int_distribution<U64> mDist;
  uniform_int_distribution<int> mSel;

  template<int N=1> U64 R64()
  {
    U64 r = mDist(mMt);
    for (int i = 1; i < N; i++)
      r &= mDist(mMt);
    return r;
  }

  U64 Magic(const int bits);
  int Transform(U64 board, const U64 magic);
  int Collisions(const U64 magic);

  static vector<U64> sBlock;
  static vector<U64> sAttack;
};

#endif // ISLAND_H

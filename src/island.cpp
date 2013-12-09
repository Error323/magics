#include "island.h"
#include <iostream>
#include <cassert>

vector<U64> Island::sAttack;
vector<U64> Island::sBlock;

void Island::Init(const int maxbits, const bool bishop, const int square)
{
  if (sAttack.empty() || sBlock.empty())
  {
    sAttack.resize(1 << maxbits);
    sBlock.resize(1 << maxbits);
    U64 mask = bishop ? bmask(square) : rmask(square);
    for (int i = 0; i < (1 << maxbits); i++)
    {
      sBlock[i] = Index2U64(i, maxbits, mask);
      sAttack[i] = bishop ? batt(square, sBlock[i]) : ratt(square, sBlock[i]);
    }
  }
}

Island::Island(int tbits, int nparents):
  mTargetBits(tbits)
{
  random_device rd;
  mMt.seed(rd());
  mSel = uniform_int_distribution<int>(0, nparents);
  mUsed.resize(1 << tbits);
}

Island::~Island()
{
}

chromosome Island::Populate()
{
  U64 best = Magic(mTargetBits);
  int best_fitness = Collisions(best);

  int fitness;
  U64 individual;
  for (int i = 1; i < num_individuals; i++)
  {
    individual = Magic(mTargetBits);
    fitness = Collisions(individual);
    if (fitness < best_fitness)
    {
      best_fitness = fitness;
      best = individual;
    }
  }

  return make_pair(best_fitness, best);
}

chromosome Island::RandChromosome()
{
  U64 magic = Magic(mTargetBits);
  return make_pair(Collisions(magic), magic);
}

chromosome Island::NextGen(const vector<chromosome> &parents, U64 best, int best_fitness)
{
  int fitness;
  U64 side, child;
  for (int i = 1; i < num_individuals; i++)
  {
    side   = R64();
    child  = parents[mSel(mMt)].second &  side;
    child |= parents[mSel(mMt)].second & ~side;
    child ^= R64<3>() & C64(0x3ffffffffffffff);

    fitness = Collisions(child);
    if (fitness <= best_fitness)
    {
      best_fitness = fitness;
      best = child;
    }
  }

  return make_pair(best_fitness, best);
}

int Island::Collisions(const U64 magic)
{
  int c = 0;
  mUsed.assign(mUsed.size(), C64(0));
  int index;

  for (int i = 0, n = sAttack.size(); i < n; i++)
  {
    index = Transform(sBlock[i], magic);
    if (mUsed[index] == C64(0))
      mUsed[index] = sAttack[i];
    else if (mUsed[index] != sAttack[i])
      c++;
  }

  return c;
}

inline U64 Island::Magic(const int bits)
{
  U64 shift = 64 - bits;
  U64 magic = R64<3>() & C64(0x3ffffffffffffff);
  magic |= shift << 58;

  return magic;
}

inline int Island::Transform(U64 board, const U64 magic)
{
  board *= magic;
  board >>= (magic >> 58);

  return static_cast<int>(board);
}

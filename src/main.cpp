#include <vector>

#include <stdlib.h>
#include <stdio.h>
#include <signal.h>
#include <math.h>
#include <getopt.h>
#include <assert.h>
#include <string.h>
#include <x86intrin.h>

#include "bitboard.h"
#include "island.h"
#include "timer.h"

constexpr int num_islands = 64;
constexpr char unit[] = {'K', 'M', 'G', 'T'};

// globals
int square;
int is_bishop;
unsigned int target_bits;
unsigned int max_bits;
unsigned int min_bits;

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

  U64 mask = is_bishop ? bmask(square) : rmask(square);
  max_bits = CountBits(mask);
  min_bits = is_bishop ? b_dist_attack_sets[square] : r_dist_attack_sets[square];
  min_bits = ceil(log2(min_bits));

  if (target_bits == 0) target_bits = max_bits;
  if (max_bits_min_one) target_bits = max_bits-1;

  if (target_bits < min_bits || target_bits > max_bits)
  {
    fprintf(stderr, "Error: target bits %d not in {%d,...,%d}\n", target_bits, min_bits, max_bits);
    return EXIT_FAILURE;
  }

  Island::Init(max_bits, is_bishop, square);

  printf("Generating magic for '%s' on square %c%d using %d <= (%d) <= %d bits\n",
          is_bishop ? "bishop" : "rook", char(square%8+65), square/8+1,
          min_bits, target_bits, max_bits);

  std::vector<Island*> islands(num_islands, nullptr);
  std::vector<chromosome> parents(num_islands+1);

  bool found = false;
  int fitness = std::numeric_limits<int>::max();
  U64 solution = C64(0);
  for (int i = 0; i < num_islands; i++)
  {
    islands[i] = new Island(target_bits, num_islands);
    parents[i] = islands[i]->Populate();
    if (parents[i].first < fitness)
    {
      solution = parents[i].second;
      fitness = parents[i].first;
      found = fitness == 0;
    }
  }

  double start_time = timer::GetRealTime();
  U32 generation = 0;
  int counter = 0;
  while (!found && !stopped)
  {
    parents.back() = islands.front()->RandChromosome();

    for (int i = 0; i < num_islands && !found; i++)
    {
      parents[i] = islands[i]->NextGen(parents, solution, fitness);
      for (int j = 0; j < i; j++)
      {
        if (parents[i].second == parents[j].second)
        {
          parents[i] = islands[i]->RandChromosome();
          break;
        }
      }

      if (parents[i].first < fitness)
      {
        solution = parents[i].second;
        fitness = parents[i].first;
        found = fitness == 0;
      }
    }

    counter++;
    double time = timer::GetRealTime() - start_time;
    if (time >= 5 || found)
    {
      for (int i = 0; i < num_islands; i++)
        fprintf(stdout, "0x%llxull %d\n", parents[i].second, parents[i].first);
      double mps = (Island::num_individuals * num_islands * counter) / time;
      int u = floor(log10(mps)) / 3;
      u = max(min(u, 4), 1);
      mps /= pow(10, u*3);
      printf("G %lu\tS %0.2f%c C %d\n", generation, mps, unit[u-1], fitness);
      start_time += time;
      counter = 0;
      fflush(stdout);
    }
    generation++;
  }

  if (found)
  {
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

  for (int i = 0; i < num_islands; i++)
    if (islands[i])
      delete islands[i];

  return EXIT_SUCCESS;
}

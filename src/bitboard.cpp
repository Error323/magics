#include "bitboard.h"

U64 Index2U64(int index, int bits, U64 m)
{
  int i, j;
  U64 result = C64(0);

  for (i = 0; i < bits; i++)
  {
    j = PopLSB(m);

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

void Print3(const U64 a, const U64 b, const U64 c)
{
  unsigned char line;
  for (int i = 7; i >= 0; i--)
  {
    printf("%d", i+1);

    line = (a >> (i * 8)) & 0xff;
    for (int j = 0; j < 8; j++)
      printf(" %s", ((line >> j) & 1) == 1 ? "1" : ".");

    printf("\t%d", i+1);

    line = (b >> (i * 8)) & 0xff;
    for (int j = 0; j < 8; j++)
      printf(" %s", ((line >> j) & 1) == 1 ? "1" : ".");

    printf("\t%d", i+1);

    line = (c >> (i * 8)) & 0xff;
    for (int j = 0; j < 8; j++)
      printf(" %s", ((line >> j) & 1) == 1 ? "1" : ".");

    printf("\n");
  }

  printf("  a b c d e f g h \t  a b c d e f g h \t  a b c d e f g h\n");
}

void Print2(const U64 a, const U64 b)
{
  unsigned char line;
  for (int i = 7; i >= 0; i--)
  {
    printf("%d", i+1);

    line = (a >> (i * 8)) & 0xff;
    for (int j = 0; j < 8; j++)
      printf(" %s", ((line >> j) & 1) == 1 ? "1" : ".");

    printf("\t%d", i+1);

    line = (b >> (i * 8)) & 0xff;
    for (int j = 0; j < 8; j++)
      printf(" %s", ((line >> j) & 1) == 1 ? "1" : ".");

    printf("\n");
  }

  printf("  a b c d e f g h \t  a b c d e f g h \n");
}

void Print1(const U64 board)
{
  for (int i = 7; i >= 0; i--)
  {
    printf("%d", i+1);
    unsigned char line = (board >> (i * 8)) & 0xff;

    for (int j = 0; j < 8; j++)
      printf(" %s", ((line >> j) & 1) == 1 ? "1" : ".");

    printf("\n");
  }

  printf("  a b c d e f g h\n");
}

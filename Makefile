CC=g++
CFLAGS= -O3 -march=core2 -Wall -Wextra

all:
	$(CC) $(CFLAGS) ga.cpp -o magics

profile:
	$(CC) $(CFLAGS) -pg ga.cpp -o magics

clean:
	git clean -f

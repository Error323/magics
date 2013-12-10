CC=g++
CFLAGS= -O3 -Wall -Wextra -std=c++11 -march=native

all:
	$(CC) $(CFLAGS) ga.cc -o magics -lrt

profile:
	$(CC) $(CFLAGS) -Wl,--no-as-needed ga.cc -o magics -lrt -lprofiler

clean:
	git clean -f

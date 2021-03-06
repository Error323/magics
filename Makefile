CC=clang++
CFLAGS= -O3 -Wall -Wextra -std=c++11 -march=native

all:
	$(CC) -o magics $(CFLAGS) ga.cc -lrt

profile:
	$(CC) -o magics $(CFLAGS) -Wl,--no-as-needed ga.cc -lrt -lprofiler

clean:
	git clean -f

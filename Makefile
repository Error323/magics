CC=g++
CFLAGS= -O2 -Wall -Wextra

all:
	$(CC) $(CFLAGS) ga.cpp -o magics

profile:
	$(CC) $(CFLAGS) -pg -g -fno-omit-frame-pointer -fno-inline-functions -fno-optimize-sibling-calls ga.cpp -o magics

debug:
	$(CC) $(CFLAGS) -g ga.cpp -o magics

clean:
	git clean -f

dot:
	gprof magics > profile.dat
	gprof2dot.py -n0 -e0 -s profile.dat > profile.dot
	xdot profile.dot

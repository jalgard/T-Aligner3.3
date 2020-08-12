CFLAGS = -Wall -O3 -std=c++1z -I./core/ -I./shell/ -Wno-unused-function -lpthread
BINDIR = ./bin
SHELLDIR = ./shell

GCC=g++

all:
	        $(GCC) $(CFLAGS) $(SHELLDIR)/alignlib.cpp -o $(BINDIR)/alignlib
	        $(GCC) $(CFLAGS) $(SHELLDIR)/mrna2ref.cpp -o $(BINDIR)/mrna2ref

	        $(GCC) $(CFLAGS) $(SHELLDIR)/snp2mrna.cpp -o $(BINDIR)/snp2mrna

test:
	        $(GCC) $(CFLAGS) $(SHELLDIR)/test.cpp -o $(BINDIR)/test

orfs:
	        $(GCC) $(CFLAGS) $(SHELLDIR)/findorfs.cpp -o $(BINDIR)/findorfs
grna:
	        $(GCC) $(CFLAGS) $(SHELLDIR)/gfinder.cpp -o $(BINDIR)/gfinder
	        $(GCC) $(CFLAGS) $(SHELLDIR)/findgrna.cpp -o $(BINDIR)/find_grna
gverify:
	        $(GCC) $(CFLAGS) $(SHELLDIR)/gverify.cpp -o $(BINDIR)/gverify

clean:
	        rm -rf *.o

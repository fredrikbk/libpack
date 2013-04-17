#!/usr/bin/make

CXX=mpic++
CC=mpicxx
MPIRUN_CMD=mpirun -n 2
FARCDIR=..

LDLIBS=-lfarcinterposer $(shell llvm-config --libs all)
LDFLAGS=-L$(FARCDIR) $(shell llvm-config --ldflags) #-dynamic
CPPFLAGS= -O3 $(shell llvm-config --cppflags)


TESTS=$(shell ls *.c | sed -e 's/.c$$//g';)

all: $(TESTS)

run: all
	@for i in $(TESTS); do $(MPIRUN_CMD) ./$$i; done

-include Makefile.deps

$(FARCDIR)/libfarc.a: $(FARCDIR)/*.cpp $(FARCDIR)/*.hpp
	make -C $(FARCDIR) libfarc.a


clean:
	rm -f $(TESTS)
	rm -f *.o
	rm -f Makefile.deps

Makefile.deps:
	$(CXX) $(CPPFLAGS) -MM *.cpp > Makefile.deps

show_ir:
	clang -emit-llvm -c input.c -o tmp.bc
	llc-3.1 -march=cpp tmp.bc -o output.cpp






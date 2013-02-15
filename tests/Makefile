CXX=clang++
CC=clang++

LDLIBS=$(shell llvm-config --libs all)
LDFLAGS=$(shell llvm-config --ldflags) -lmpich
CPPFLAGS=-g3 $(shell llvm-config --cppflags) -I/usr/include/mpi 

BINARIES=$(shell ls *.cpp | sed -e 's/.cpp//';)

all: $(BINARIES)
	./runtests.sh

-include Makefile.deps

clean:
	rm -f $(BINARIES)
	rm -f *.o
	rm -f Makefile.deps

Makefile.deps:
	$(CXX) $(CPPFLAGS) -MM *.cpp > Makefile.deps

show_ir:
	clang -emit-llvm -c input.c -o tmp.bc
	llc-3.1 -march=cpp tmp.bc -o output.cpp






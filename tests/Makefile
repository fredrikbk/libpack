CXX=mpic++
CC=mpicc

LDLIBS=$(shell llvm-config --libs all)
LDFLAGS=$(shell llvm-config --ldflags) #-dynamic
CPPFLAGS= -O3 $(shell llvm-config --cppflags) -I../..

BINARIES=$(shell ls *.cpp | sed -e 's/.cpp$$//g';)

all: $(BINARIES)
	./runtests.sh

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CPPFLAGS)

%: %.o ../ddt_jit.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

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






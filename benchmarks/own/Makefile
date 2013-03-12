CXX=mpic++

LDLIBS=$(shell llvm-config --libs all)
LDFLAGS=$(shell llvm-config --ldflags) -lmpich -dynamic
CPPFLAGS=-DHRT_ARCH=2 -O3 $(shell llvm-config --cppflags) -I/usr/include/mpi 

BINARIES=$(shell ls *.cpp | sed -e 's/.cpp//';)

all: $(BINARIES)

clean:
	rm -f $(BINARIES)

show_ir:
	clang -emit-llvm -c input.c -o tmp.bc
	llc-3.1 -march=cpp tmp.bc -o output.cpp

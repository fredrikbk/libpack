CXX=clang++

LDLIBS=$(shell llvm-config --libs all)
LDFLAGS=$(shell llvm-config --ldflags) -lmpich
CPPFLAGS=-g3 $(shell llvm-config --cppflags) -I/usr/include/mpi 

all: interposer.a

-include Makefile.deps

Makefile.deps:
	$(CXX) $(CPPFLAGS) -MM *.cpp > Makefile.deps

interposer.a: interposer.o
	 ar r $@ interposer.o
	 ranlib $@

clean:
	rm -f *.o interposer.a Makefile.deps

show_ir:
	clang -emit-llvm -c input.c -o tmp.bc
	llc -march=cpp tmp.bc -o output.cpp


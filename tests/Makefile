CXX=mpic++
CC=mpicc

FARCDIR=..

LDLIBS=-lfarc $(shell llvm-config --libs all)
LDFLAGS=$(shell llvm-config --ldflags) -L$(FARCDIR) #-dynamic
CPPFLAGS= -O3 $(shell llvm-config --cppflags) -I../..


TESTS = $(shell ls *.cpp | sed -e 's/.cpp$$//g';)
TESTS += $(shell ls *.c | sed -e 's/.c$$//g';)

all: $(TESTS)

run: all
	@for i in $(TESTS); do ./$$i; done

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CPPFLAGS)

%.o: %.c
	$(CC) -c -o $@ $< $(CPPFLAGS)

%: %.o $(FARCDIR)/libfarc.a
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)


# delete the default rules which compile everything in one go
%: %.cpp

%: %.c

$(FARCDIR)/libfarc.a: $(FARCDIR)/*.cpp $(FARCDIR)/*.hpp
	make -C $(FARCDIR) libfarc.a

-include Makefile.deps

clean:
	rm -f $(TESTS)
	rm -f *.o
	rm -f Makefile.deps

Makefile.deps:
	$(CXX) $(CPPFLAGS) -MM *.cpp > Makefile.deps

show_ir:
	clang -emit-llvm -c input.c -o tmp.bc
	llc-3.1 -march=cpp tmp.bc -o output.cpp






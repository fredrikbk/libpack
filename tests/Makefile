CXX=mpic++
CC=mpicc

LDLIBS=$(shell llvm-config --libs all)
LDFLAGS=$(shell llvm-config --ldflags) #-dynamic
CPPFLAGS= -O3 $(shell llvm-config --cppflags) -I../..

FARCDIR=..

TESTS=$(shell ls *.cpp | sed -e 's/.cpp$$//g';)

all: $(TESTS)

run: all
	@for i in $(TESTS); do ./$$i; done

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CPPFLAGS)

%: %.o $(FARCDIR)/libfarc.a
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

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






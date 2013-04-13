MPICXX=mpic++
CXX=clang++
LEX=flex

.PHONY: all

FARCDIR=..

LDLIBS=$(shell llvm-config --libs all)
LDFLAGS=$(shell llvm-config --ldflags) #-dynamic
CPPFLAGS=-DHRT_ARCH=2 -O3 $(shell llvm-config --cppflags) -I/usr/include/mpi -I$(FARCDIR) -I$(FARCDIR)/copy_benchmark/hrtimer

OBJECTS=lexer.o parser.o $(FARCDIR)/libfarc.a

ddtplayer: $(OBJECTS)
	$(MPICXX) -o $@ $^ $(LDFLAGS) $(LDLIBS) 

$(FARCDIR)/libfarc.a: $(FARCDIR)/*.cpp $(FARCDIR)/*.hpp
	make -C $(FARCDIR) libfarc.a

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CPPFLAGS)

parser.cpp parser.hpp: parser.y
	bison -d -o parser.cpp $<

lexer.c: lexer.flex parser.cpp
	flex -o $@ $<

clean:
	rm -f ddtplayer *.o parser.hpp parser.cpp lexer.c


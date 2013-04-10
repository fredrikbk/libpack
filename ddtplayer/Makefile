MPICXX=mpic++
CXX=clang++
LEX=flex

.PHONY: all parser

FARCDIR=..

LDLIBS=$(shell llvm-config --libs all)
LDFLAGS=$(shell llvm-config --ldflags) #-dynamic
CPPFLAGS=-DHRT_ARCH=2 -O3 $(shell llvm-config --cppflags) -I/usr/include/mpi -I$(FARCDIR) -I$(FARCDIR)/copy_benchmark/hrtimer

OBJECTS=lexer.o parser.o $(FARCDIR)/libfarc.a

ddtplayer: $(OBJECTS)
	$(MPICXX) -o $@ $^ $(LDFLAGS) $(LDLIBS) 

all: parser ddtplayer

$(FARCDIR)/libfarc.a: $(FARCDIR)/*.cpp $(FARCDIR)/*.hpp
	make -C $(FARCDIR) libfarc.a

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CPPFLAGS)

parser: parser.y lexer.lex
	bison -d -o parser.cpp parser.y
	flex -o lexer.c lexer.lex

clean:
	rm -f ddtplayer *.o


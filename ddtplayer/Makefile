MPICXX=mpic++
CXX=clang++
LEX=flex

.PHONY: all

FARCDIR=..

LDLIBS=$(shell llvm-config --libs all)
LDFLAGS=$(shell llvm-config --ldflags) #-dynamic
CPPFLAGS=-DHRT_ARCH=2 -O3 $(shell llvm-config --cppflags) -I/usr/include/mpi -I$(FARCDIR) -I$(FARCDIR)/copy_benchmark/hrtimer

OBJECTS=cmdline.o ddtplayer.o lex.yy.o parser.tab.o cmdline.o $(FARCDIR)/libfarc.a

ddtplayer: $(OBJECTS) 
	$(MPICXX) -o $@ $^ $(LDFLAGS) $(LDLIBS) 

$(FARCDIR)/libfarc.a: $(FARCDIR)/*.cpp $(FARCDIR)/*.hpp
	make -C $(FARCDIR) libfarc.a

parser.tab.cpp parser.tab.hpp: parser.ypp
	bison -d $<

lex.yy.c: lexer.flex parser.tab.cpp
	flex $<

cmdline.h cmdline.c: cmdline.ggo
	gengetopt < cmdline.ggo

clean:
	rm -f ddtplayer *.o  *.tab.cpp *.tab.hpp lex.yy.c cmdline.h cmdline.c


#!/usr/bin/make

CXX=clang++

CFLAGS_MPI=-D_FORTIFY_SOURCE=2 -I/usr/include/mpich2
LDFLAGS_MPI=-Wl,-z,relro -L/usr/lib
LDLIBS_MPI=-lmpich -lopa -lrt -lpthread # -lmpl -lcr

LDLIBS=$(shell llvm-config --libs all) $(LDLIBS_MPI)
LDFLAGS=$(shell llvm-config --ldflags) $(LDFLAGS_MPI)
CFLAGS=-g3 $(shell llvm-config --cppflags) $(CFLAGS_MPI)

BINARIES=$(shell ls *.c | sed -e 's/.c$$//g';)

SOURCES=$(shell ls *.c;) ../interposer.cpp ../interposer_common.cpp

all: $(BINARIES)
	./runtests.sh

-include Makefile.deps

Makefile.deps:
	$(CXX) $(CFLAGS) -MM -x c++ $(SOURCES) > Makefile.deps

interposer.o:
	$(CXX) $(CFLAGS) -c ../interposer.cpp -o interposer.o

interposer_common.o:
	$(CXX) $(CFLAGS) -c ../interposer_common.cpp -o interposer_common.o

%.o: %.c
	$(CC) -c $(CFLAGS) $< -o $@

%: %.o interposer.o interposer_common.o
	$(CXX) $(LDFLAGS) $< interposer.o interposer_common.o -o $@ $(LDLIBS)

%: %.c

clean:
	rm -f *.o
	rm -f Makefile.deps
	rm -f $(BINARIES)


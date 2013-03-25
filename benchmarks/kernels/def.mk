CXX=g++
FARCDIR=../../..

FARCLIBS=$(shell llvm-config --libs all)
FARCLDFLAGS=$(shell llvm-config --ldflags)

CPPFLAGS= -O3 -march=core-avx-i -I/usr/include/mpi -I$(FARCDIR) $(shell llvm-config --cppflags)

.SUFFIXES:
.PRECIOUS: %.s %.o %.ll

%.o: %.s
	$(CXX) $(CPPFLAGS) -c -o $@ $< 

%.s: %.cpp
	$(CXX) $(CPPFLAGS) -S -o $@ $< 

%.ll: %.cpp
	clang $(CPPFLAGS) -S -emit-llvm -O3 -o $@ $< 

%: %.o ../driver.o
	$(CXX) -o $@ $^ 

$(FARCDIR)/ddt_jit.o: $(FARCDIR)/ddt_jit.cpp $(FARCDIR)/ddt_jit.hpp
	@make -C $(FARCDIR) ddt_jit.o

%_farc: %_farc.o ../driver.o $(FARCDIR)/ddt_jit.o
	$(CXX) -o $@ $^ $(FARCLDFLAGS) $(FARCLIBS)

%_mpi: %_mpi.o ../driver.o
	mpic++ -o $@ $^ 

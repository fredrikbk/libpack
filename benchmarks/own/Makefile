CXX=mpic++

LDLIBS=$(shell llvm-config --libs all)
LDFLAGS=$(shell llvm-config --ldflags) #-dynamic
CPPFLAGS=-DHRT_ARCH=2 -O3 -g3 $(shell llvm-config --cppflags) -I/usr/include/mpi -I../.. 

FARC=-L../.. -lfarcinterposer

BINARIES=$(shell ls *.cpp | sed -e 's/.cpp//';)

all: $(BINARIES)

../../ddt_jit.o: ../../ddt_jit.cpp ../../ddt_jit.hpp
	make -C ../../ ddt_jit.o

%: %.o ../../ddt_jit.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CPPFLAGS)

vector_pack: vector_pack.o ../../ddt_jit.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

speedtest_vector_pack: speedtest_vector_pack.o ../../ddt_jit.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

speedtest_vector_send: speedtest_vector_send.o ../../interposer.o ../../interposer_common.o ../../ddt_jit.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

speedtest_indexed_block_pack: speedtest_indexed_block_pack.o ../../ddt_jit.o
	$(CXX) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(BINARIES) *.o

show_ir:
	clang -emit-llvm -c input.c -o tmp.bc
	llc-3.1 -march=cpp tmp.bc -o output.cpp

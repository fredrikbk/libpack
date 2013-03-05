CXX=mpic++
F77=mpif77

LDLIBS=$(shell llvm-config --libs all)
LDFLAGS=$(shell llvm-config --ldflags)
CFLAGS=$(shell llvm-config --cppflags)

all: lib
	make -C tests
	make -C pmpi-tests

lib: interposer.o interposer_common.o
	ar rcs libfarc.a interposer.o interposer_common.o 

flib: interposerf interposer_common.o
	ar rcs libfarcf.a interposerf.o interposer_common.o 

interposer.o: interposer.c
	$(CXX) $(CFLAGS) -c interposer.c -o interposer.o

interposer_common.o: interposer_common.cpp
	$(CXX) $(CFLAGS) -c interposer_common.cpp -o interposer_common.o

interposerf: interposer.f
	$(F77) -c interposer.f -o interposerf.o

clean:
	rm -f *.o *.a
	make -C tests clean
	make -C pmpi-tests clean
	svn status

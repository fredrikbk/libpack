CXX=mpic++
F77=mpif77

LDLIBS=$(shell llvm-config --libs all)
LDFLAGS=$(shell llvm-config --ldflags)
CFLAGS=$(shell llvm-config --cppflags)

farc: interposer.o interposerf.o interposer_common.o

test: farc
	make -C tests
	make -C pmpi-tests

interposer.o: interposer.c
	$(CXX) -c interposer.c -o interposer.o

interposerf.o: interposerf.f
	$(F77) -c interposerf.f -o interposerf.o

interposer_common.o: interposer_common.cpp
	$(CXX) $(CFLAGS) -DHRT_ARCH=2 -c interposer_common.cpp -o interposer_common.o

clean:
	rm -f *.o
	make -C tests clean
	make -C pmpi-tests clean
	svn status

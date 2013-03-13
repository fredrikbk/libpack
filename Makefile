CXX=mpic++
F77=mpif77

LDLIBS=$(shell llvm-config --libs all)
LDFLAGS=$(shell llvm-config --ldflags)
CFLAGS=$(shell llvm-config --cppflags)

#farc: libfarc.a libfarcinterposer.a

farc: ddt_jit.o interposer.o interposerf.o interposer_common.o

test: farc
	make -C tests
	make -C pmpi-tests

#libfarc.a: ddt_jit.o
#	ar rcs libfarc.a $^

#libfarcinterposer.a: interposer.o interposer_common.o ddt_jit.o
#	ar rcs libfarcinterposer.a $^

interposer.o: interposer.c interposer_common.h
	$(CXX) -DHRT_ARCH=2 -c $< -o $@

interposerf.o: interposerf.f interposer_common.h
	$(F77) -c $< -o $@

interposer_common.o: interposer_common.cpp ddt_jit.hpp
	$(CXX) $(CFLAGS) -DHRT_ARCH=2 -c $< -o $@

ddt_jit.o: ddt_jit.cpp
	$(CXX) $(CFLAGS) -DHRT_ARCH=2 -c ddt_jit.cpp -o ddt_jit.o


clean:
	rm -f *.o libfarcinterposer.a
	make -C tests clean
	make -C pmpi-tests clean
	svn status

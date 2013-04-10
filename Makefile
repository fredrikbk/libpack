CXX=mpic++
F77=mpif77

# Configuration variables picked up from the environment
PACKVAR     ?= 0
LLVM_OUTPUT ?= 0

CONFIGVARS = -DPACKVAR=$(PACKVAR) -DLLVM_OUTPUT=$(LLVM_OUTPUT)

FARC= ddt_jit.o codegen_common.o codegen_primitive.o codegen_contiguous.o codegen_vector.o codegen_indexed.o

LDLIBS+=$(shell llvm-config --libs all)
LDFLAGS+=$(shell llvm-config --ldflags)
CPPFLAGS+=$(shell llvm-config --cppflags) -O3 -g3

farc: libfarcinterposer.a libfarc.a

test: farc
	make run -C tests
	make -C pmpi-tests

libfarc.a: $(FARC)
	ar rcs libfarc.a $^

libfarcinterposer.a: interposer.o interposer_common.o $(FARC)
	ar rcs libfarcinterposer.a $^

interposer.o: interposer.c interposer_common.h
	$(CXX) -DHRT_ARCH=2 -c $< -o $@

interposerf.o: interposerf.f interposer_common.h
	$(F77) -c $< -o $@

interposer_common.o: interposer_common.cpp ddt_jit.hpp
	$(CXX) $(CPPFLAGS) -DHRT_ARCH=2 -c $< -o $@

%.o: %.cpp codegen.hpp codegen_common.hpp ddt_jit.hpp
	$(CXX) -c -o $@ $< $(CPPFLAGS) $(CONFIGVARS) -DHRT_ARCH=2


clean:
	rm -f *.o *.a
	make -C tests clean
	make -C pmpi-tests clean
	svn status

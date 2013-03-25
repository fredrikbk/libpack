FC=mpif90
CC=mpicc
LD=mpicxx

PAPI_INSTALL_PATH=
FARC_PATH=/home/timos/Work/farc-mpiddt

HRT_ARCH=6   # 1 : x86_32, gcc, rdtsc
             # 2 : x64_64, gcc, rdtsc
             # 3 : ppc, gcc, tb
             # 4 : ia64, gcc, itc
             # 5 : mips64-sicortex, gcc
             # 6 : mpi_wtime

TEST_TYPE=1  # 1: time
             # 2: papi
             # 3: time+papi

OPTIMIZATION=-O3
DEBUG=
WARNINGS=-Wall

FCFLAGS+=${OPTIMIZATION} ${DEBUG} ${WARNINGS} -DTEST_TYPE=${TEST_TYPE} -g
CCFLAGS+=${OPTIMIZATION} ${DEBUG} ${WARNINGS} -DHRT_ARCH=${HRT_ARCH} -DTEST_TYPE=${TEST_TYPE} -g
LDFLAGS=$(shell llvm-config --ldflags) -L${FARC_PATH} 
LIBS=-lpapi -lfarcinterposer $(shell llvm-config --libs all)
F90_LIBS=-lgfortran

ifneq ($(PAPI_INSTALL_PATH),)
	CCFLAGS+=-I$(PAPI_INSTALL_PATH)/include
	LDFLAGS+=-L$(PAPI_INSTALL_PATH)/lib
endif

FC=ftn
CC=cc
LD=CC

PAPI_INSTALL_PATH=
FARC_PATH=/u/staff/timos/Work/farc-mpiddt

HRT_ARCH=2   # 1 : x86_32, gcc, rdtsc
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

FCFLAGS+=${OPTIMIZATION} ${DEBUG} ${WARNINGS} -DTEST_TYPE=${TEST_TYPE} -dynamic
CCFLAGS+=${OPTIMIZATION} ${DEBUG} ${WARNINGS} -DHRT_ARCH=${HRT_ARCH} -DTEST_TYPE=${TEST_TYPE} -dynamic
LDFLAGS=$(shell llvm-config --ldflags) -L${FARC_PATH} -dynamic

# without the farc wrapper
LIBS=-lpapi $(shell llvm-config --libs all)

# #with the farc wrapper
#LIBS=-lpapi -lfarcinterposer $(shell llvm-config --libs all)

F90_LIBS=-lgfortran

ifneq ($(PAPI_INSTALL_PATH),)
	CCFLAGS+=-I$(PAPI_INSTALL_PATH)/include
	LDFLAGS+=-L$(PAPI_INSTALL_PATH)/lib
endif

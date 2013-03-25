#!/bin/bash

execdir=`pwd`

#get the latest papi
if [ ! -d "papi-5.0.1/src" ]; then
  wget 'http://icl.cs.utk.edu/projects/papi/downloads/papi-5.0.1.tar.gz'
  tar xzvf papi-5.0.1.tar.gz
  cd papi-5.0.1/src
  ./configure --prefix=$execdir/papi_inst
  make
  make install
  cd ../..
fi
sed -i "s+PAPI_INSTALL_PATH=.*$+PAPI_INSTALL_PATH=$execdir/papi_inst+" Makefile.inc
export LD_LIBRARY_PATH="$execdir/papi_inst/lib:$LD_LIBRARY_PATH"

#sed -i "s/HRT_ARCH=2/HRT_ARCH=6/" Makefile.inc

# meassure only time in the first run

sed -i "s/TEST_TYPE=[0-9]/TEST_TYPE=1/" Makefile.inc

make clean
make

salloc -N 2 --ntasks-per-node=1 mpirun -n 2 ./ddtbench_c
mv ddtbench.out ddtbench_c_time_inter.out

salloc -N 2 --ntasks-per-node=1 mpirun -n 2 ./ddtbench_f90
mv ddtbench.out ddtbench_f90_time_inter.out

salloc -N 1 --ntasks-per-node=2 mpirun -n 2 ./ddtbench_c
mv ddtbench.out ddtbench_c_time_intra.out

salloc -N 1 --ntasks-per-node=2 mpirun -n 2 ./ddtbench_f90
mv ddtbench.out ddtbench_f90_time_intra.out

# Meassure time and PAPI counters

# sed -i "s/TEST_TYPE=[0-9]/TEST_TYPE=3/" Makefile.inc
# 
# make clean
# make
# 
# export PAPI_EVT1="PAPI_L1_TCM"
# export PAPI_EVT2="PAPI_L2_TCM"
# 
# mpirun -n 2 ./ddtbench_c
# mv ddtbench.out ddtbench_c_papi.out
# 
# mpirun -n 2 ./ddtbench_f90
# mv ddtbench.out ddtbench_f90_papi.out




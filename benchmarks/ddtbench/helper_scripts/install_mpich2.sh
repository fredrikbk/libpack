#!/bin/sh

install=$1

mkdir -p local/$install
cd local/$install
TEMP=`pwd`
cd ../../src_mpi/
if [ ! -d "$install" ]; then
  tar xvfz $install.tar.gz
fi
cd $install
./configure --prefix=$TEMP
make
make install
cd ../..

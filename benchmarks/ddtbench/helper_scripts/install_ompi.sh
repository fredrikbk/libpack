#!/bin/sh

install=$1

mkdir -p local/$install
cd local/$install
TEMP=`pwd`
cd ../../src_mpi/
if [ ! -d "$install" ]; then
  tar xvfj $install.tar.bz2
fi
cd $install
./configure --prefix=$TEMP
make
make install
cd ../..

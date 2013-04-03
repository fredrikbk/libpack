#!/bin/bash
set -e

if [ -z "$1" ]
then
    echo "gencode.sh <output-name> <count> <blocklen> <stride>"
    exit 1
fi

if [ -z "$2" ]
then
    echo "gencode.sh <output-name> <count> <blocklen> <stride>"
    exit 1
fi

if [ -z "$3" ]
then
    echo "gencode.sh <output-name> <count> <blocklen> <stride>"
    exit 1
fi

if [ -z "$4" ]
then
    echo "gencode.sh <output-name> <count> <blocklen> <stride>"
    exit 1
fi

rm -f ../../ddt_jit.o
LLVM_OUTPUT=1 make vector_pack
rm -f ../../ddt_jit.o

N="$1-$2-$3-$4.ll"

./vector_pack $2 $3 $4 2> $N
llc -O3 $N

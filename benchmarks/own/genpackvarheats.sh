#!/bin/bash
for i in `seq 1 8`
do
    rm -f ../../ddt_jit.o
    PACKVAR=$i ./genheat.sh "pv$i"
done

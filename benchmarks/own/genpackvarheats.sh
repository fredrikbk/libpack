#!/bin/bash
for i in `seq 1 9`
do
    rm -f ../../ddt_jit.o
    PACKVAR=$i ./genheat.sh "pv$i"
done

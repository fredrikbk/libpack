#!/bin/bash
for i in `seq 1 9`
do
    rm ../../ddt_jit.o
    PACKVAR=$i ./genheat_light.sh "pv$i"
done

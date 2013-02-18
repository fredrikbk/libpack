#!/bin/sh

TESTS=`ls test_* | grep -v '\.'`

for i in $TESTS; do
    mpirun -np 2 ./$i 2>/dev/null
done

exit 0

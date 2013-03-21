#!/bin/sh

TESTS=`ls test_*`

for i in $TESTS; do
    ./$i
done

exit 0

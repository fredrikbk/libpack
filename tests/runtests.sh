#!/bin/sh

TESTS=`ls test_*`

for i in $TESTS; do
    ./$i 2>/dev/null
done

exit 0

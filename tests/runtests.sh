#!/bin/sh

TESTS=`find . -name "test_*" -executable | sort`

for i in $TESTS; do
    ./$i
    # 2>/dev/null
done

exit 0

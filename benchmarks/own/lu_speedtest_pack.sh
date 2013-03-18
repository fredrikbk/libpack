#!/bin/bash
pushd ../../; make; popd
make speedtest_vector_pack && ./speedtest_vector_pack 5 5 5 1 180 180 1 32 32 1 1 1 1

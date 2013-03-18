#!/bin/bash
make speedtest_vector_send && mpirun -n 2 ./speedtest_vector_send 5 5 5 1 180 180 1 32 32 1 1 1 1

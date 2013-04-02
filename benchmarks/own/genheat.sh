#!/bin/bash
set -e

if [ -z "$1" ]
then
    echo "genheat.sh <output-name>"
    exit 1
fi

make speedtest_vector_pack
./speedtest_vector_pack 30 8 80 8 180 180 1 10 1000 20 1 1 1 > speedtest_vector_pack.dat
cat heatmap.R | R --no-save

mv speedtest_vector_pack.dat "$1.dat"
mv speedtest_vector_pack.svg "$1.svg"


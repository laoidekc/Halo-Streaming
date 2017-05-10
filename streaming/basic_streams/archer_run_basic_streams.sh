#!/bin/bash
make clean
bolt -n $1 -N $2 -t 0:5:0 -o basic_streams.bolt -j basic_streams -A d35-clk ./basic_streams 100000 2 1000000 $1
make
qsub basic_streams.bolt
watch qstat -u $USER

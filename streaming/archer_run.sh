#!/bin/bash
make clean
bolt -n $1 -N $2 -t 0:5:0 -o halo_streaming.bolt -j halo_streaming -A d68 ./halo_streaming
make archer
qsub halo_streaming.bolt $1
watch qstat -u $USER
